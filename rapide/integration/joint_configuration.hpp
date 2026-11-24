#pragma once
#include "rapide/integration/configuration.hpp"
#include "rapide/model/joint.hpp"

namespace rapide::integration {

/** Advance packed joint coordinates; free-joint velocity is body-expressed. */
template <typename S, typename Q, typename V>
inline Eigen::Matrix<S, Eigen::Dynamic, 1>
integrate_joint_configuration(const JointCoordinates<S> &joint,
                              const Eigen::MatrixBase<Q> &configuration,
                              const Eigen::MatrixBase<V> &velocity, S dt) {
  if (configuration.size() != joint_nq(joint) ||
      velocity.size() != joint_nv(joint)) {
    throw std::invalid_argument("Joint configuration/velocity size mismatch");
  }
  Eigen::Matrix<S, Eigen::Dynamic, 1> next = configuration;
  if (joint_type(joint) == JointType::Free) {
    const Eigen::Quaternion<S> orientation(configuration(3), configuration(0),
                                           configuration(1), configuration(2));
    const auto next_orientation =
        integrate_orientation(orientation, velocity.template head<3>(), dt);
    const Eigen::Matrix<S, 3, 1> world_velocity =
        orientation * velocity.template tail<3>();
    next.template head<4>() = next_orientation.coeffs();
    next.template tail<3>() =
        advance_linear(configuration.template tail<3>(), world_velocity, dt);
  } else if (joint_nv(joint) != 0) {
    next(0) = advance_linear(configuration(0), velocity(0), dt);
  }
  return next;
}

} // namespace rapide::integration
