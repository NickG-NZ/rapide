#pragma once
#include "rapide/math/rotation.hpp"

namespace rapide::integration {

/** Advance a Euclidean quantity with a derivative held constant over the step.
 */
template <typename Value, typename Rate, typename S>
inline typename Value::PlainObject
advance_linear(const Eigen::MatrixBase<Value> &value,
               const Eigen::MatrixBase<Rate> &rate, S dt) {
  return (value + rate * dt).eval();
}

template <typename S> inline S advance_linear(S value, S rate, S dt) {
  return value + rate * dt;
}

/** Right multiplication applies a delta expressed in the body's current axes.
 */
template <typename S>
inline Eigen::Quaternion<S>
apply_body_rotation(const Eigen::Quaternion<S> &orientation,
                    const Eigen::Quaternion<S> &delta) {
  return (orientation * delta).normalized();
}

/** Exponential orientation update, exact for constant body angular velocity. */
template <typename S, typename Derived>
inline Eigen::Quaternion<S>
integrate_orientation(const Eigen::Quaternion<S> &orientation,
                      const Eigen::MatrixBase<Derived> &body_angular_velocity,
                      S dt) {
  return apply_body_rotation(orientation,
                             rotation_increment(body_angular_velocity * dt));
}

} // namespace rapide::integration
