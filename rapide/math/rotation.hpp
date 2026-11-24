#pragma once
#include "rapide/math/linalg.hpp"
#include <cmath>

namespace rapide {

template <typename Derived>
inline Eigen::Quaternion<typename Derived::Scalar>
rotation_increment(const Eigen::MatrixBase<Derived> &angle) {
  using S = typename Derived::Scalar;
  using std::cos;
  using std::sin;
  using std::sqrt;
  const S squared_angle = angle.squaredNorm();
  if (squared_angle < S(1e-20)) {
    return Eigen::Quaternion<S>(S(1), angle.x() / S(2), angle.y() / S(2),
                                angle.z() / S(2))
        .normalized();
  }
  const S n = sqrt(squared_angle);
  const S scale = sin(n / S(2)) / n;
  return Eigen::Quaternion<S>(cos(n / S(2)), scale * angle.x(),
                              scale * angle.y(), scale * angle.z());
}

inline Vec3 rotation_log(Eigen::Quaterniond q) {
  q.normalize();
  if (q.w() < 0) {
    q.coeffs() *= -1;
  }
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

} // namespace rapide
