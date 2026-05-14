#pragma once
#include "rapide/model/model.hpp"
#include "rapide/model/joint.hpp"
#include "rapide/model/body.hpp"
#include "rapide/math/se3.hpp"

namespace rapide {
namespace vehicles {

// Motorcycle model parameters (SI units).
template<typename S>
struct MotorcycleParams {
    // Chassis
    S chassis_mass   = S(180);
    Eigen::Matrix<S,3,1> chassis_com{S(0), S(0), S(0.5)};   // above chassis ref point
    Eigen::Matrix<S,3,3> chassis_I = (Eigen::Matrix<S,3,3>() <<
        S(20),  S(0), S(0),
        S(0),  S(80), S(0),
        S(0),   S(0), S(80)
    ).finished();

    // Wheels (identical front/rear in v1)
    S wheel_mass     = S(15);
    S wheel_radius   = S(0.3);
    Eigen::Matrix<S,3,3> wheel_I = (Eigen::Matrix<S,3,3>() <<
        S(0.5), S(0), S(0),
        S(0),  S(1),  S(0),
        S(0),   S(0), S(1)
    ).finished();

    // Geometry
    S wheelbase      = S(1.4);   // rear-to-front wheel centre distance
    S chassis_height = S(0.5);   // chassis reference point height above wheel axes
};

// build_motorcycle<S>() — construct kinematic tree for v1 motorcycle.
//
// Topology:
//   world → chassis [FreeJoint]
//           ├── rear_wheel  [RevoluteJoint, axis +Y]
//           └── front_wheel [RevoluteJoint, axis +Y]
//
// The chassis reference frame origin is at the rear wheel axle.
template<typename S>
Model<S> build_motorcycle(const MotorcycleParams<S>& p = {}) {
    ModelBuilder<S> builder;

    // ── Chassis (floating base) ───────────────────────────────────────────────
    // Attached to world via FreeJoint. Initial position: chassis at height h
    // above ground so wheels touch at z=0.
    Body<S> chassis = Body<S>::with_inertia(
        "chassis", p.chassis_mass, p.chassis_com, p.chassis_I);

    // Fixed transform from world to joint predecessor = identity (joint is at origin).
    int chassis_idx = builder.add_body(
        std::move(chassis),
        FreeJoint<S>{},
        /*parent=*/0,
        SE3<S>::Identity());

    // ── Rear wheel ───────────────────────────────────────────────────────────
    // Joint predecessor frame: at the rear axle position in chassis frame.
    // Chassis ref origin IS at the rear axle, so transform is identity offset upward
    // to chassis CoM height ... actually the chassis ref origin is at the rear axle.
    // The rear axle is at (0, 0, 0) in chassis frame.
    Body<S> rear_wheel = Body<S>::with_inertia(
        "rear_wheel", p.wheel_mass,
        Eigen::Matrix<S,3,1>::Zero(),  // CoM at wheel centre
        p.wheel_I);

    Eigen::Matrix<S,3,1> rear_axle_in_chassis(S(0), S(0), S(0));
    builder.add_body(
        std::move(rear_wheel),
        RevoluteJoint<S>{Eigen::Matrix<S,3,1>(S(0), S(1), S(0))},  // +Y axis
        chassis_idx,
        SE3<S>{Eigen::Matrix<S,3,3>::Identity(), rear_axle_in_chassis});

    // ── Front wheel ──────────────────────────────────────────────────────────
    // Front axle is at (wheelbase, 0, 0) in chassis frame (no fork in v1).
    Body<S> front_wheel = Body<S>::with_inertia(
        "front_wheel", p.wheel_mass,
        Eigen::Matrix<S,3,1>::Zero(),
        p.wheel_I);

    Eigen::Matrix<S,3,1> front_axle_in_chassis(p.wheelbase, S(0), S(0));
    builder.add_body(
        std::move(front_wheel),
        RevoluteJoint<S>{Eigen::Matrix<S,3,1>(S(0), S(1), S(0))},
        chassis_idx,
        SE3<S>{Eigen::Matrix<S,3,3>::Identity(), front_axle_in_chassis});

    return builder.build();
}

// Convenience: build a default initial state with the motorcycle standing upright.
// Chassis centre at height = wheel_radius above ground.
template<typename S>
State<S> motorcycle_standing_state(const Model<S>& model, S wheel_radius = S(0.3)) {
    State<S> s = model.zero_state();
    // FreeJoint q slice: [qx, qy, qz, qw, px, py, pz]
    // Identity rotation (qw=1) already set by zero_state().
    // Set z height so wheels touch the ground.
    s.q(6) = wheel_radius;  // pz = wheel_radius (chassis ref at rear axle height)
    return s;
}

}  // namespace vehicles
}  // namespace rapide
