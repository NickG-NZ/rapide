#pragma once
#include <Eigen/Core>
#include <Eigen/Geometry>

namespace rapide {

using Vec3 = Eigen::Vector3d;
using Mat3 = Eigen::Matrix3d;

/** Cross-product matrix: skew(v) * w == v.cross(w). */
template <typename Derived>
inline Eigen::Matrix<typename Derived::Scalar, 3, 3>
skew(const Eigen::MatrixBase<Derived> &v) {
  EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived, 3);
  using S = typename Derived::Scalar;
  Eigen::Matrix<S, 3, 3> result;
  result << S(0), -v(2), v(1), v(2), S(0), -v(0), -v(1), v(0), S(0);
  return result;
}

} // namespace rapide
