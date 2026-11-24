#pragma once
#include "rapide/contact/tire.hpp"
#include "rapide/simulator.hpp"
#include <numbers>

namespace rapide::examples {

template <typename S> struct MotorcycleParams {
  // Chassis
  S chassis_mass = S(180);
  Eigen::Matrix<S, 3, 3> chassis_I =
      (Eigen::Matrix<S, 3, 3>() << S(20), S(0), S(0), S(0), S(80), S(0), S(0),
       S(0), S(80))
          .finished();
  Eigen::Matrix<S, 3, 1> chassis_com{S(0.7), S(0),
                                     S(0.1)}; // forward of rear axle, above

  // Wheels (identical front/rear for v1)
  S wheel_mass = S(15);
  S wheel_radius = S(0.3);
  Eigen::Matrix<S, 3, 3> wheel_I =
      (Eigen::Matrix<S, 3, 3>() << S(0.5), S(0), S(0), S(0), S(1.0), S(0), S(0),
       S(0), S(1.0))
          .finished();

  // Geometry
  S wheelbase = S(1.4); // rear-to-front axle distance

  // Suspension (only used when has_suspension = true)
  bool has_suspension = false;
  S fork_mass = S(4);        // mass of each slider/fork body [kg]
  S suspension_k = S(25000); // spring stiffness [N/m]
  S suspension_c = S(1500);  // damping [N·s/m]

  // Tire contact
  S tire_k_n = S(1e5); // normal stiffness [N/m]
  S tire_c_n = S(1e3); // normal damping [N·s/m]
  S tire_mu = S(0.8);  // friction coefficient
};

template <typename S> struct MotorcycleModel {
  Model model;
  BodyId rear_wheel_body, front_wheel_body;
  JointId rear_wheel_joint, front_wheel_joint;
  JointId rear_fork_joint = -1, front_fork_joint = -1;
};

template <typename S>
MotorcycleModel<S> build_motorcycle(const MotorcycleParams<S> &p = {}) {
  ModelBuilder builder;
  const double radius = double(p.wheel_radius);
  const auto chassis = builder.add_body(
      "chassis", double(p.chassis_mass), p.chassis_I.template cast<double>(),
      {{0, 0, radius}}, p.chassis_com.template cast<double>());
  const auto hinge = rotation_increment(Vec3{-std::numbers::pi / 2, 0, 0});
  MotorcycleModel<S> result;
  for (int i = 0; i < 2; ++i) {
    const double x = i ? double(p.wheelbase) : 0;
    BodyId parent = chassis;
    Pose wheel_frame{{x, 0, 0}, hinge};
    JointId fork = -1;
    if (p.has_suspension) {
      parent = builder.add_body(i ? "front_slider" : "rear_slider",
                                double(p.fork_mass), .01 * Mat3::Identity(),
                                {{x, 0, radius}});
      fork = builder.add_joint(JointType::Prismatic, chassis, parent,
                               {{x, 0, 0}}, {});
      wheel_frame.position.setZero();
    }
    const auto wheel =
        builder.add_body(i ? "front_wheel" : "rear_wheel", double(p.wheel_mass),
                         p.wheel_I.template cast<double>(), {{x, 0, radius}});
    const auto joint = builder.add_joint(JointType::Revolute, parent, wheel,
                                         wheel_frame, {{0, 0, 0}, hinge});
    if (i) {
      result.front_wheel_body = wheel;
      result.front_wheel_joint = joint;
      result.front_fork_joint = fork;
    } else {
      result.rear_wheel_body = wheel;
      result.rear_wheel_joint = joint;
      result.rear_fork_joint = fork;
    }
  }
  result.model = builder.build();
  return result;
}

template <typename S>
MotorcycleModel<S> build_motorcycle_suspended(MotorcycleParams<S> p = {}) {
  p.has_suspension = true;
  return build_motorcycle(p);
}

inline State motorcycle_standing_state(const Model &model,
                                       double wheel_radius = .3) {
  auto state = model.initial_state;
  const auto &chassis = state.bodies.at(1);
  const double origin_height =
      (chassis.pose.position - chassis.pose.orientation * model.bodies[1].com)
          .z();
  for (size_t i = 0; i < model.bodies.size(); ++i) {
    if (model.bodies[i].dynamic()) {
      state.bodies[i].pose.position.z() += wheel_radius - origin_height;
    }
  }
  return state;
}

template <typename S>
void attach_tire_forces(Simulator &sim, const MotorcycleModel<S> &bike,
                        const MotorcycleParams<S> &p) {
  const TireContact<double> tire{double(p.wheel_radius), double(p.tire_k_n),
                                 double(p.tire_c_n), double(p.tire_mu)};
  sim.set_force_callback([=](const Model &model, const State &state) {
    auto forces = zero_forces(model);
    for (auto body : {bike.rear_wheel_body, bike.front_wheel_body}) {
      const auto &s = state.bodies[body];
      const auto force =
          tire.compute({s.pose.orientation.toRotationMatrix(), s.pose.position},
                       {s.omega, s.pose.orientation.conjugate() * s.velocity});
      forces[body] = {force.torque, force.force};
    }
    return forces;
  });
}

template <typename S>
void attach_suspension(Simulator &sim, const MotorcycleModel<S> &bike,
                       const MotorcycleParams<S> &p) {
  for (auto joint : {bike.rear_fork_joint, bike.front_fork_joint}) {
    if (joint >= 0) {
      sim.add_joint_spring(joint, double(p.suspension_k),
                           double(p.suspension_c), 0);
    }
  }
}

} // namespace rapide::examples
