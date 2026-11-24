#pragma once
#include "rapide/model/body.hpp"
#include "rapide/model/joint.hpp"
#include "rapide/model/state.hpp"
#include <Eigen/Cholesky>

namespace rapide {
using JointId = int;
struct Model {
  std::vector<Body> bodies;
  std::vector<Joint> joints;
  State initial_state;
  int nv =
      0; // Size of the shared packed body-velocity vector, not ABA joint DOFs.

  void validate() const {
    if (bodies.empty() || bodies[0].dynamic()) {
      throw std::invalid_argument(
          "Model needs a static world body at index zero");
    }
    int velocity_size = 0;
    for (const auto &body : bodies) {
      if (!body.com.allFinite()) {
        throw std::invalid_argument("Nonfinite COM");
      }
      if (!body.dynamic()) {
        continue;
      }
      Eigen::LLT<Mat3> inertia(body.inertia);
      if (!std::isfinite(body.mass) || body.mass <= 0 ||
          !body.inertia.allFinite() ||
          !body.inertia.isApprox(body.inertia.transpose(), 1e-12) ||
          inertia.info() != Eigen::Success ||
          body.velocity_offset != velocity_size ||
          !body.inverse_inertia.allFinite() ||
          !(body.inertia * body.inverse_inertia)
               .isApprox(Mat3::Identity(), 1e-8)) {
        throw std::invalid_argument(
            "Invalid body mass, inertia or velocity layout");
      }
      velocity_size += 6;
    }
    if (nv != velocity_size) {
      throw std::invalid_argument("Invalid body velocity size");
    }
    for (const auto &joint : joints) {
      joint.validate();
      if (joint.a >= static_cast<int>(bodies.size()) ||
          joint.b >= static_cast<int>(bodies.size()) ||
          (!bodies[joint.a].dynamic() && !bodies[joint.b].dynamic())) {
        throw std::invalid_argument("Invalid joint body pair");
      }
    }
    validate_state(initial_state);
  }

  void validate_state(const State &state) const {
    if (state.bodies.size() != bodies.size()) {
      throw std::invalid_argument("State does not match model");
    }
    for (size_t i = 0; i < bodies.size(); ++i) {
      const auto &s = state.bodies[i];
      validate_pose(s.pose);
      if (!s.omega.allFinite() || !s.velocity.allFinite()) {
        throw std::invalid_argument("Velocity must be finite");
      }
      if (!bodies[i].dynamic()) {
        const auto &fixed = initial_state.bodies[i];
        if ((s.pose.position - fixed.pose.position).norm() > 1e-12 ||
            s.pose.orientation.angularDistance(fixed.pose.orientation) >
                1e-12 ||
            s.omega.norm() != 0 || s.velocity.norm() != 0) {
          throw std::invalid_argument("Static body state cannot change");
        }
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
        llt.info() != Eigen::Success || !com.allFinite()) {
      throw std::invalid_argument(
          "Dynamic mass and COM inertia must be finite and positive definite");
    }
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
        b >= static_cast<int>(model_.bodies.size()) || a == b) {
      throw std::invalid_argument("Joint needs two distinct valid body IDs");
    }
    if (!model_.bodies[a].dynamic() && !model_.bodies[b].dynamic()) {
      throw std::invalid_argument("Joint needs at least one dynamic body");
    }
    const Joint joint{type, a, b, frame_a, frame_b};
    joint.validate();
    if (type == JointType::Free) {
      throw std::invalid_argument("Free motion needs no joint");
    }
    model_.joints.push_back(joint);
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

} // namespace rapide
