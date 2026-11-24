#pragma once
#include "rapide/math/pose.hpp"
#include <string>

namespace rapide {
using BodyId = int;
struct Body {
  std::string name;
  double mass = 0;
  Mat3 inertia = Mat3::Zero();
  Mat3 inverse_inertia = Mat3::Zero();
  Vec3 com = Vec3::Zero();
  int velocity_offset = -1;
  bool dynamic() const { return velocity_offset >= 0; }
};

} // namespace rapide
