#pragma once
#include "rapide/math/pose.hpp"
#include <vector>

namespace rapide {
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

} // namespace rapide
