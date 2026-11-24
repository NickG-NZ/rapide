#pragma once
#include "rapide/math/rotation.hpp"
#include "rapide/model/model.hpp"

namespace rapide {

/** World-expressed force and torque about the body's COM. */
struct BodyForce {
  Vec3 torque = Vec3::Zero();
  Vec3 force = Vec3::Zero();
};

using Forces = std::vector<BodyForce>;

struct LinearSpring {
  BodyId a, b;
  Vec3 point_a, point_b;
  double stiffness, damping, rest_length;
};

inline Forces zero_forces(const Model &model) {
  return Forces(model.bodies.size());
}

inline void add_force_at_point(const Model &model, const State &state,
                               Forces &forces, BodyId body, const Vec3 &local,
                               const Vec3 &force) {
  auto &f = forces.at(body);
  const Vec3 arm = state.bodies.at(body).pose.orientation *
                   (local - model.bodies.at(body).com);
  f.force += force;
  f.torque += arm.cross(force);
}

/** Apply equal-and-opposite spring/damper forces between authored local points.
 */
inline void add_spring(const Model &model, const State &state, Forces &forces,
                       BodyId a, const Vec3 &point_a, BodyId b,
                       const Vec3 &point_b, double stiffness, double damping,
                       double rest_length) {
  if (!std::isfinite(stiffness) || !std::isfinite(damping) ||
      !std::isfinite(rest_length) || stiffness < 0 || damping < 0 ||
      rest_length < 0) {
    throw std::invalid_argument("Invalid spring parameters");
  }
  const Vec3 d = world_point(model, state, b, point_b) -
                 world_point(model, state, a, point_a);
  const double length = d.norm();
  if (length < 1e-12) {
    throw std::runtime_error(
        "Coincident spring endpoints have no force direction");
  }
  const Vec3 n = d / length;
  const double rate = n.dot(point_velocity(model, state, b, point_b) -
                            point_velocity(model, state, a, point_a));
  const Vec3 force = (stiffness * (length - rest_length) + damping * rate) * n;
  add_force_at_point(model, state, forces, a, point_a, force);
  add_force_at_point(model, state, forces, b, point_b, -force);
}

inline void add_spring(const Model &model, const State &state, Forces &forces,
                       const LinearSpring &spring) {
  add_spring(model, state, forces, spring.a, spring.point_a, spring.b,
             spring.point_b, spring.stiffness, spring.damping,
             spring.rest_length);
}

/** Revolute angle is wrapped to [-pi,pi]; prismatic coordinate is signed
 * travel. */
inline double joint_position(const Model &model, const State &state,
                             JointId id) {
  const auto &j = model.joints.at(id);
  const Eigen::Quaterniond qa =
      state.bodies[j.a].pose.orientation * j.frame_a.orientation;
  if (j.type == JointType::Revolute) {
    const Eigen::Quaterniond qb =
        state.bodies[j.b].pose.orientation * j.frame_b.orientation;
    return rotation_log(qa.conjugate() * qb).z();
  }
  if (j.type == JointType::Prismatic) {
    return (qa * Vec3::UnitZ())
        .dot(world_point(model, state, j.b, j.frame_b.position) -
             world_point(model, state, j.a, j.frame_a.position));
  }
  throw std::invalid_argument("Joint has no scalar coordinate");
}

inline double joint_velocity(const Model &model, const State &state,
                             JointId id) {
  const auto &j = model.joints.at(id);
  const auto &a = state.bodies[j.a];
  const auto &b = state.bodies[j.b];
  const Vec3 axis = a.pose.orientation * j.frame_a.orientation * Vec3::UnitZ();
  const Vec3 wa = a.pose.orientation * a.omega;
  if (j.type == JointType::Revolute) {
    return axis.dot(b.pose.orientation * b.omega - wa);
  }
  if (j.type == JointType::Prismatic) {
    const Vec3 d = world_point(model, state, j.b, j.frame_b.position) -
                   world_point(model, state, j.a, j.frame_a.position);
    return axis.dot(point_velocity(model, state, j.b, j.frame_b.position) -
                    point_velocity(model, state, j.a, j.frame_a.position)) +
           wa.cross(axis).dot(d);
  }
  throw std::invalid_argument("Joint has no scalar velocity");
}

/** Positive effort acts on body B along/about frame A's +Z axis. */
inline void add_joint_effort(const Model &model, const State &state,
                             Forces &forces, JointId id, double effort) {
  if (!std::isfinite(effort)) {
    throw std::invalid_argument("Nonfinite joint effort");
  }
  const auto &j = model.joints.at(id);
  const Vec3 axis = state.bodies[j.a].pose.orientation * j.frame_a.orientation *
                    Vec3::UnitZ();
  if (j.type == JointType::Revolute) {
    forces.at(j.a).torque -= effort * axis;
    forces.at(j.b).torque += effort * axis;
  } else if (j.type == JointType::Prismatic) {
    add_force_at_point(model, state, forces, j.a, j.frame_a.position,
                       -effort * axis);
    add_force_at_point(model, state, forces, j.b, j.frame_b.position,
                       effort * axis);
  } else {
    throw std::invalid_argument("Joint does not accept scalar effort");
  }
}

inline void add_joint_spring(const Model &model, const State &state,
                             Forces &forces, JointId id, double stiffness,
                             double damping, double rest) {
  if (!std::isfinite(stiffness) || !std::isfinite(damping) ||
      !std::isfinite(rest) || stiffness < 0 || damping < 0) {
    throw std::invalid_argument("Invalid joint spring parameters");
  }
  add_joint_effort(model, state, forces, id,
                   -stiffness * (joint_position(model, state, id) - rest) -
                       damping * joint_velocity(model, state, id));
}

} // namespace rapide
