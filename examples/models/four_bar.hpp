#pragma once
#include "rapide/math/pose.hpp"
#include "rapide/model/forces.hpp"

namespace rapide::examples {

struct FourBar {
  rapide::Model model;
  rapide::BodyId crank, coupler, follower;
  rapide::JointId drive;
  rapide::BodyId chassis;
};

/** Planar parallelogram represented by four spatial hinges, including a loop.
 *  Three moving bodies have 18 velocities and 20 rows of rank 17. */
inline FourBar build_four_bar(double angle = .6,
                              const rapide::Pose &mount = {}) {
  ModelBuilder builder;
  const BodyId chassis = builder.add_static("chassis", mount);
  const auto mounted = [&](const Pose &pose) -> Pose {
    return {mount.position + mount.orientation * pose.position,
            mount.orientation * pose.orientation};
  };
  const Vec3 direction{std::cos(angle), std::sin(angle), 0};
  const auto rotation = rotation_increment(Vec3{0, 0, angle});
  const Mat3 crank_inertia = Vec3{.01, 1. / 12, 1. / 12}.asDiagonal();
  const BodyId a = builder.add_body("crank", 1, crank_inertia,
                                    mounted({.5 * direction, rotation}));
  const BodyId b =
      builder.add_body("coupler", 2, Vec3{.02, 2. / 3, 2. / 3}.asDiagonal(),
                       mounted({direction + Vec3{1, 0, 0}}));
  const BodyId c =
      builder.add_body("follower", 1, crank_inertia,
                       mounted({Vec3{2, 0, 0} + .5 * direction, rotation}));
  const JointId drive =
      builder.add_joint(JointType::Revolute, chassis, a, {}, {{-.5, 0, 0}});
  builder.add_joint(JointType::Revolute, a, b, {{.5, 0, 0}}, {{-1, 0, 0}});
  builder.add_joint(JointType::Revolute, b, c, {{1, 0, 0}}, {{.5, 0, 0}});
  builder.add_joint(JointType::Revolute, chassis, c, {{2, 0, 0}},
                    {{-.5, 0, 0}});
  return {builder.build(), a, b, c, drive, chassis};
}

/** Unloaded linear shock between the chassis and coupler at the initial pose.
 */
inline rapide::LinearSpring four_bar_spring(const FourBar &mechanism) {
  LinearSpring spring{
      mechanism.chassis, mechanism.coupler, {1.5, 0, 0}, {0, 0, 0}, 500, 50, 0};
  const auto &m = mechanism.model;
  spring.rest_length =
      (world_point(m, m.initial_state, spring.b, spring.point_b) -
       world_point(m, m.initial_state, spring.a, spring.point_a))
          .norm();
  return spring;
}

} // namespace rapide::examples
