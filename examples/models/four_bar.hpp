#pragma once
#include "rapide/maximal/model.hpp"

namespace rapide::examples {

struct FourBar {
  maximal::Model model;
  maximal::BodyId crank, coupler, follower;
  maximal::JointId drive;
};

/** Planar parallelogram represented by four spatial hinges, including a loop.
 *  Three moving bodies have 18 velocities and 20 rows of rank 17. */
inline FourBar build_four_bar(double angle = .6) {
  using namespace maximal;
  ModelBuilder builder;
  const Vec3 direction{std::cos(angle), std::sin(angle), 0};
  const auto rotation = rotation_increment(Vec3{0, 0, angle});
  const Mat3 crank_inertia = Vec3{.01, 1. / 12, 1. / 12}.asDiagonal();
  const BodyId a =
      builder.add_body("crank", 1, crank_inertia, {.5 * direction, rotation});
  const BodyId b =
      builder.add_body("coupler", 2, Vec3{.02, 2. / 3, 2. / 3}.asDiagonal(),
                       {direction + Vec3{1, 0, 0}});
  const BodyId c = builder.add_body("follower", 1, crank_inertia,
                                    {Vec3{2, 0, 0} + .5 * direction, rotation});
  const JointId drive =
      builder.add_joint(JointType::Revolute, 0, a, {}, {{-.5, 0, 0}});
  builder.add_joint(JointType::Revolute, a, b, {{.5, 0, 0}}, {{-1, 0, 0}});
  builder.add_joint(JointType::Revolute, b, c, {{1, 0, 0}}, {{.5, 0, 0}});
  builder.add_joint(JointType::Revolute, 0, c, {{2, 0, 0}}, {{-.5, 0, 0}});
  return {builder.build(), a, b, c, drive};
}

} // namespace rapide::examples
