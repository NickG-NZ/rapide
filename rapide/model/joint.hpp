#pragma once
#include "rapide/math/pose.hpp"
#include "rapide/math/se3.hpp"
#include "rapide/math/spatial.hpp"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <variant>

namespace rapide {

enum class JointType { Ball, Fixed, Revolute, Prismatic, Free };

/** Physical attachment frames relative to authored body origins; axis is +Z. */
struct Joint {
  JointType type;
  int a, b;
  Pose frame_a, frame_b;

  void validate() const {
    switch (type) {
    case JointType::Ball:
    case JointType::Fixed:
    case JointType::Revolute:
    case JointType::Prismatic:
    case JointType::Free:
      break;
    default:
      throw std::invalid_argument("Unknown joint type");
    }
    if (a < 0 || b < 0 || a == b) {
      throw std::invalid_argument("Joint needs two distinct valid body IDs");
    }
    validate_pose(frame_a);
    validate_pose(frame_b);
  }
};

// ─── RevoluteJoint ───────────────────────────────────────────────────────────
// Single-axis rotation. nq = nv = 1.
template <typename S> struct RevoluteJoint {
  static constexpr JointType type = JointType::Revolute;
  Eigen::Matrix<S, 3, 1> axis; // unit rotation axis in predecessor frame

  static constexpr int NQ = 1;
  static constexpr int NV = 1;

  SE3<S> X_J(S q) const {
    Eigen::AngleAxis<S> aa(q, axis);
    return {aa.toRotationMatrix(), Eigen::Matrix<S, 3, 1>::Zero()};
  }

  // Motion subspace: S_col = [axis; 0]
  Twist<S> S_times_v(S dq) const {
    return {axis * dq, Eigen::Matrix<S, 3, 1>::Zero()};
  }
  Twist<S> vJ(S dq) const { return S_times_v(dq); }

  // S^T * f = axis · f.torque
  S S_transpose_force(const Wrench<S> &f) const { return axis.dot(f.torque); }
};

// ─── PrismaticJoint ──────────────────────────────────────────────────────────
// Single-axis translation. nq = nv = 1.
template <typename S> struct PrismaticJoint {
  static constexpr JointType type = JointType::Prismatic;
  Eigen::Matrix<S, 3, 1> axis; // unit translation axis in predecessor frame

  static constexpr int NQ = 1;
  static constexpr int NV = 1;

  SE3<S> X_J(S q) const {
    return {Eigen::Matrix<S, 3, 3>::Identity(), axis * q};
  }

  // Motion subspace: S_col = [0; axis]
  Twist<S> S_times_v(S dq) const {
    return {Eigen::Matrix<S, 3, 1>::Zero(), axis * dq};
  }
  Twist<S> vJ(S dq) const { return S_times_v(dq); }

  // S^T * f = axis · f.force   (translational — NOT torque)
  S S_transpose_force(const Wrench<S> &f) const { return axis.dot(f.force); }
};

// ─── FreeJoint
// ──────────────────────────────────────────────────────────────── 6-DOF
// floating base. nq = 7 (quaternion xyzw + position xyz), nv = 6. Velocity is
// body-frame spatial velocity [ω; v].
template <typename S> struct FreeJoint {
  static constexpr JointType type = JointType::Free;
  static constexpr int NQ = 7;
  static constexpr int NV = 6;

  SE3<S> X_J(const Eigen::Ref<const Eigen::Matrix<S, 7, 1>> &q) const {
    Eigen::Quaternion<S> quat(q(3), q(0), q(1), q(2)); // (w, x, y, z)
    quat.normalize();
    return {quat.toRotationMatrix(), q.template tail<3>()};
  }

  Twist<S> vJ(const Eigen::Ref<const Eigen::Matrix<S, 6, 1>> &dq) const {
    return {dq.template head<3>(), dq.template tail<3>()};
  }

  Eigen::Matrix<S, 6, 1> S_transpose_force(const Wrench<S> &f) const {
    Eigen::Matrix<S, 6, 1> tau;
    tau << f.torque, f.force;
    return tau;
  }
};

template <typename S> struct FixedJoint {
  static constexpr JointType type = JointType::Fixed;
  static constexpr int NQ = 0;
  static constexpr int NV = 0;
  SE3<S> X_J() const { return SE3<S>::Identity(); }
};

// ─── ABA joint coordinates ───────────────────────────────────────────────────
template <typename S>
using JointCoordinates = std::variant<FreeJoint<S>, RevoluteJoint<S>,
                                      PrismaticJoint<S>, FixedJoint<S>>;

template <typename S> JointType joint_type(const JointCoordinates<S> &joint) {
  return std::visit([](const auto &value) { return value.type; }, joint);
}

template <typename S> int joint_nq(const JointCoordinates<S> &j) {
  return std::visit([](const auto &jt) { return jt.NQ; }, j);
}

template <typename S> int joint_nv(const JointCoordinates<S> &j) {
  return std::visit([](const auto &jt) { return jt.NV; }, j);
}

} // namespace rapide
