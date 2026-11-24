#pragma once
#include "rapide/math/linalg.hpp"
#include <cmath>
#include <stdexcept>

namespace rapide {

/** Rigid pose: position and orientation expressed in a parent frame. */
struct Pose {
  Vec3 position = Vec3::Zero();
  Eigen::Quaterniond orientation = Eigen::Quaterniond::Identity();
};

inline void validate_pose(const Pose &p) {
  if (!p.position.allFinite() || !p.orientation.coeffs().allFinite() ||
      std::abs(p.orientation.norm() - 1) > 1e-8) {
    throw std::invalid_argument("Pose must be finite with a unit quaternion");
  }
}

} // namespace rapide
