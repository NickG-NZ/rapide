#pragma once
#include "rapide/contact/compliant_tire.hpp"
#include <numbers>

namespace rapide::examples {

struct Bike {
  rapide::Model model;
  rapide::BodyId chassis, rear, front, swingarm;
  rapide::JointId suspension;
  rapide::LinearSpring spring;
  rapide::Tire tire;
};

/** Single-pivot rear rocker with an opposing spring arm; steering is fixed. */
inline Bike build_bike(double forward_speed = 2, double lateral_speed = 0,
                       double drop_height = 0) {
  if (!std::isfinite(drop_height) || drop_height < 0) {
    throw std::invalid_argument("Drop height must be finite and nonnegative");
  }
  Tire tire;
  tire.stiffness = 10000;
  tire.damping = 200;
  tire.longitudinal_slip_stiffness = 300;
  tire.lateral_slip_stiffness = 300;
  constexpr double frame_mass = 30, wheel_mass = 3, arm_mass = 2;
  const double total_mass = frame_mass + 2 * wheel_mass + arm_mass;
  const double rear_load = 9.81 * (total_mass * .65 + arm_mass * .2) / 1.3;
  const double front_load = total_mass * 9.81 - rear_load;
  const double height = tire.radius - rear_load / tire.stiffness;
  const double front_height = tire.radius - front_load / tire.stiffness;
  ModelBuilder builder;
  const auto chassis =
      builder.add_body("bike_frame", frame_mass, Vec3{1, 4, 4}.asDiagonal(),
                       {{0, 0, height + .32}});
  const double inertia = .5 * wheel_mass * tire.radius * tire.radius;
  const Mat3 wheel_inertia =
      Vec3{inertia / 2, inertia, inertia / 2}.asDiagonal();
  const auto rear = builder.add_body("rear_wheel", wheel_mass, wheel_inertia,
                                     {{-.65, 0, height}});
  const auto front = builder.add_body("front_wheel", wheel_mass, wheel_inertia,
                                      {{.65, 0, front_height}});
  const auto hinge = rotation_increment(Vec3{-std::numbers::pi / 2, 0, 0});
  const auto swingarm = builder.add_body("rear_swingarm", arm_mass,
                                         Vec3{.03, .15, .15}.asDiagonal(),
                                         {{-.2, 0, height + .15}});
  const auto suspension =
      builder.add_joint(JointType::Revolute, chassis, swingarm,
                        {{-.2, 0, -.17}, hinge}, {{0, 0, 0}, hinge});
  builder.add_joint(JointType::Revolute, swingarm, rear,
                    {{-.45, 0, -.15}, hinge}, {{0, 0, 0}, hinge});
  builder.add_joint(JointType::Revolute, chassis, front,
                    {{.65, 0, front_height - height - .32}, hinge},
                    {{0, 0, 0}, hinge});
  LinearSpring spring{chassis, swingarm, {.05, 0, -.42}, {.25, 0, .3}, 5000,
                      150,     0};
  const double preload = .45 * (rear_load - wheel_mass * 9.81) / .25;
  spring.rest_length = .55 + preload / spring.stiffness;
  auto model = builder.build();
  for (auto body : {chassis, rear, front, swingarm}) {
    model.initial_state.bodies[body].velocity =
        Vec3{forward_speed, lateral_speed, 0};
  }
  for (auto body : {rear, front}) {
    model.initial_state.bodies[body].omega.y() =
        forward_speed / model.initial_state.bodies[body].pose.position.z();
  }
  for (auto body : {chassis, rear, front, swingarm}) {
    model.initial_state.bodies[body].pose.position.z() += drop_height;
  }
  return {std::move(model), chassis,    rear,   front,
          swingarm,         suspension, spring, tire};
}

inline rapide::Forces bike_forces(const Bike &bike, const rapide::Model &model,
                                  const rapide::State &state) {
  auto forces = rapide::zero_forces(model);
  for (auto wheel : {bike.rear, bike.front}) {
    forces[wheel] = bike.tire.compute(state.bodies[wheel]).wrench;
  }
  rapide::add_spring(model, state, forces, bike.spring);
  return forces;
}

} // namespace rapide::examples
