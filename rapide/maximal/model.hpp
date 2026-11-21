#pragma once
#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace rapide::maximal {

using Vec3 = Eigen::Vector3d;
using Mat3 = Eigen::Matrix3d;
using BodyId = int;
using JointId = int;

inline Mat3 skew(const Vec3 &x) {
  Mat3 m;
  m << 0, -x.z(), x.y(), x.z(), 0, -x.x(), -x.y(), x.x(), 0;
  return m;
}

inline Eigen::Quaterniond rotation_increment(const Vec3 &angle) {
  const double n = angle.norm();
  if (n < 1e-10)
    return Eigen::Quaterniond(1, angle.x() / 2, angle.y() / 2, angle.z() / 2)
        .normalized();
  return Eigen::Quaterniond(Eigen::AngleAxisd(n, angle / n));
}

inline Vec3 rotation_log(Eigen::Quaterniond q) {
  q.normalize();
  if (q.w() < 0)
    q.coeffs() *= -1;
  const double n = q.vec().norm();
  return n < 1e-10 ? Vec3(2 * q.vec())
                   : Vec3(2 * std::atan2(n, q.w()) / n * q.vec());
}

inline Mat3 inverse_left_jacobian(const Vec3 &angle) {
  const double n = angle.norm();
  const Mat3 k = skew(angle);
  const double c = n < 1e-5 ? 1. / 12 + n * n / 720
                            : (1 - .5 * n / std::tan(.5 * n)) / (n * n);
  return Mat3::Identity() - .5 * k + c * k * k;
}

struct Pose {
  Vec3 position = Vec3::Zero();
  Eigen::Quaterniond orientation = Eigen::Quaterniond::Identity();
};

inline void validate_pose(const Pose &p) {
  if (!p.position.allFinite() || !p.orientation.coeffs().allFinite() ||
      std::abs(p.orientation.norm() - 1) > 1e-8)
    throw std::invalid_argument("Pose must be finite with a unit quaternion");
}

/** Pose position is the COM in world coordinates; omega is body-expressed,
 *  relative to world, and linear velocity is world-expressed at the COM. */
struct BodyState {
  Pose pose;
  Vec3 omega = Vec3::Zero();
  Vec3 velocity = Vec3::Zero();
};

struct State {
  std::vector<BodyState> bodies;
};

struct Body {
  std::string name;
  double mass = 0;
  Mat3 inertia = Mat3::Zero();
  Mat3 inverse_inertia = Mat3::Zero();
  Vec3 com = Vec3::Zero();
  int velocity_offset = -1;
  bool dynamic() const { return velocity_offset >= 0; }
};

enum class JointType { Ball, Fixed, Revolute, Prismatic };

/** Frames are relative to authored body origins. The joint axis is local +Z. */
struct Joint {
  JointType type;
  BodyId a, b;
  Pose frame_a, frame_b;
  int rows() const {
    return type == JointType::Ball ? 3 : type == JointType::Fixed ? 6 : 5;
  }
};

struct Model {
  std::vector<Body> bodies;
  std::vector<Joint> joints;
  State initial_state;
  int nv = 0;

  void validate_state(const State &state) const {
    if (state.bodies.size() != bodies.size())
      throw std::invalid_argument("State does not match model");
    for (size_t i = 0; i < bodies.size(); ++i) {
      const auto &s = state.bodies[i];
      validate_pose(s.pose);
      if (!s.omega.allFinite() || !s.velocity.allFinite())
        throw std::invalid_argument("Velocity must be finite");
      if (!bodies[i].dynamic()) {
        const auto &fixed = initial_state.bodies[i];
        if ((s.pose.position - fixed.pose.position).norm() > 1e-12 ||
            s.pose.orientation.angularDistance(fixed.pose.orientation) >
                1e-12 ||
            s.omega.norm() != 0 || s.velocity.norm() != 0)
          throw std::invalid_argument("Static body state cannot change");
      }
    }
  }
};

class ModelBuilder {
public:
  ModelBuilder() { add_static("world"); }

  /** Add a dynamic body with inertia about its COM, expressed in body axes.
   *  authored_pose locates the body origin, not its COM. */
  BodyId add_body(std::string name, double mass, const Mat3 &inertia,
                  const Pose &authored_pose = {},
                  const Vec3 &com = Vec3::Zero()) {
    validate_pose(authored_pose);
    Eigen::LLT<Mat3> llt(inertia);
    if (!std::isfinite(mass) || mass <= 0 || !inertia.allFinite() ||
        !inertia.isApprox(inertia.transpose(), 1e-12) ||
        llt.info() != Eigen::Success || !com.allFinite())
      throw std::invalid_argument(
          "Dynamic mass and COM inertia must be finite and positive definite");
    Body body{std::move(name), mass, inertia, llt.solve(Mat3::Identity()), com,
              model_.nv};
    model_.nv += 6;
    model_.bodies.push_back(std::move(body));
    Pose p = authored_pose;
    p.position += p.orientation * com;
    model_.initial_state.bodies.push_back({p});
    return static_cast<int>(model_.bodies.size() - 1);
  }

  BodyId add_static(std::string name, const Pose &pose = {}) {
    validate_pose(pose);
    Body body;
    body.name = std::move(name);
    model_.bodies.push_back(std::move(body));
    model_.initial_state.bodies.push_back({pose});
    return static_cast<int>(model_.bodies.size() - 1);
  }

  JointId add_joint(JointType type, BodyId a, BodyId b,
                    const Pose &frame_a = {}, const Pose &frame_b = {}) {
    if (a < 0 || b < 0 || a >= static_cast<int>(model_.bodies.size()) ||
        b >= static_cast<int>(model_.bodies.size()) || a == b)
      throw std::invalid_argument("Joint needs two distinct valid body IDs");
    if (!model_.bodies[a].dynamic() && !model_.bodies[b].dynamic())
      throw std::invalid_argument("Joint needs at least one dynamic body");
    validate_pose(frame_a);
    validate_pose(frame_b);
    model_.joints.push_back({type, a, b, frame_a, frame_b});
    return static_cast<int>(model_.joints.size() - 1);
  }

  Model build() const { return model_; }

private:
  Model model_;
};

inline Vec3 world_point(const Model &model, const State &state, BodyId body,
                        const Vec3 &local) {
  const auto &s = state.bodies.at(body);
  return s.pose.position +
         s.pose.orientation * (local - model.bodies.at(body).com);
}

inline Vec3 point_velocity(const Model &model, const State &state, BodyId body,
                           const Vec3 &local) {
  const auto &s = state.bodies.at(body);
  return s.velocity + s.pose.orientation *
                          (s.omega.cross(local - model.bodies.at(body).com));
}

inline Eigen::VectorXd velocities(const Model &model, const State &state) {
  Eigen::VectorXd v(model.nv);
  for (size_t i = 0; i < model.bodies.size(); ++i) {
    const int offset = model.bodies[i].velocity_offset;
    if (offset < 0)
      continue;
    v.segment<3>(offset) = state.bodies[i].omega;
    v.segment<3>(offset + 3) = state.bodies[i].velocity;
  }
  return v;
}

inline void set_velocities(const Model &model, State &state,
                           const Eigen::VectorXd &v) {
  if (v.size() != model.nv || !v.allFinite())
    throw std::invalid_argument("Invalid velocity vector");
  for (size_t i = 0; i < model.bodies.size(); ++i) {
    const int offset = model.bodies[i].velocity_offset;
    if (offset < 0)
      continue;
    state.bodies[i].omega = v.segment<3>(offset);
    state.bodies[i].velocity = v.segment<3>(offset + 3);
  }
}

/** Apply body-axis rotation and world-axis translation increments. Preserve
 *  world angular velocity when correcting a body's orientation. */
inline void retract(const Model &model, State &state,
                    const Eigen::VectorXd &delta) {
  if (delta.size() != model.nv || !delta.allFinite())
    throw std::invalid_argument("Invalid pose increment");
  for (size_t i = 0; i < model.bodies.size(); ++i) {
    const int offset = model.bodies[i].velocity_offset;
    if (offset < 0)
      continue;
    auto &s = state.bodies[i];
    const Eigen::Quaterniond dq = rotation_increment(delta.segment<3>(offset));
    s.pose.orientation = (s.pose.orientation * dq).normalized();
    s.pose.position += delta.segment<3>(offset + 3);
    s.omega = dq.conjugate() * s.omega;
  }
}

} // namespace rapide::maximal
