#pragma once
#include "rapide/contact/compliant_tire.hpp"

namespace rapide::examples {

struct PlanarWheel {
  rapide::Model model;
  rapide::BodyId wheel;
  rapide::Tire tire;
};

/** Symmetric initial conditions and forces preserve XZ motion without clamps.
 */
inline PlanarWheel build_planar_wheel(double height = .47, double speed = .65,
                                      double camber = 0) {
  ModelBuilder builder;
  constexpr double mass = 12, radius = .35;
  const double inertia = .5 * mass * radius * radius;
  const auto wheel = builder.add_body(
      "wheel", mass, Vec3{inertia / 2, inertia, inertia / 2}.asDiagonal(),
      {{0, 0, height}, rotation_increment(Vec3{camber, 0, 0})});
  auto model = builder.build();
  model.initial_state.bodies[wheel].velocity.x() = speed;
  return {std::move(model), wheel, {}};
}

} // namespace rapide::examples
