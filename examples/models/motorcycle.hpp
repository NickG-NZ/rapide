#pragma once
#include "rapide/model/model.hpp"
#include "rapide/model/joint.hpp"
#include "rapide/model/body.hpp"
#include "rapide/math/se3.hpp"
#include "rapide/contact/tire.hpp"
#include "rapide/contact/suspension.hpp"
#include "rapide/simulator.hpp"

namespace rapide {
namespace examples {

// ─── Parameters ──────────────────────────────────────────────────────────────

template<typename S>
struct MotorcycleParams {
    // Chassis
    S chassis_mass = S(180);
    Eigen::Matrix<S,3,3> chassis_I = (Eigen::Matrix<S,3,3>() <<
        S(20), S(0), S(0),
        S(0), S(80), S(0),
        S(0), S(0), S(80)
    ).finished();
    Eigen::Matrix<S,3,1> chassis_com{S(0.7), S(0), S(0.1)};  // forward of rear axle, above

    // Wheels (identical front/rear for v1)
    S wheel_mass   = S(15);
    S wheel_radius = S(0.3);
    Eigen::Matrix<S,3,3> wheel_I = (Eigen::Matrix<S,3,3>() <<
        S(0.5), S(0), S(0),
        S(0), S(1.0), S(0),
        S(0), S(0), S(1.0)
    ).finished();

    // Geometry
    S wheelbase    = S(1.4);   // rear-to-front axle distance

    // Suspension (only used when has_suspension = true)
    bool has_suspension = false;
    S fork_mass        = S(4);    // mass of each slider/fork body [kg]
    S suspension_k     = S(25000); // spring stiffness [N/m]
    S suspension_c     = S(1500);  // damping [N·s/m]

    // Tire contact
    S tire_k_n = S(1e5);   // normal stiffness [N/m]
    S tire_c_n = S(1e3);   // normal damping [N·s/m]
    S tire_mu  = S(0.8);   // friction coefficient
};

// ─── Build result ─────────────────────────────────────────────────────────────

template<typename S>
struct MotorcycleModel {
    Model<S> model;
    int rear_wheel_body  = -1;
    int front_wheel_body = -1;
    int rear_fork_joint  = -1;   // PrismaticJoint index (suspension only)
    int front_fork_joint = -1;   // PrismaticJoint index (suspension only)
};

// ─── Topology without suspension ──────────────────────────────────────────────
//
//   world → chassis [FreeJoint]
//           ├── rear_wheel  [RevoluteJoint, +Y]
//           └── front_wheel [RevoluteJoint, +Y]
//
template<typename S>
MotorcycleModel<S> build_motorcycle(const MotorcycleParams<S>& p = {}) {
    ModelBuilder<S> builder;

    // Chassis
    int chassis_idx = builder.add_body(
        Body<S>::with_inertia("chassis", p.chassis_mass, p.chassis_com, p.chassis_I),
        FreeJoint<S>{},
        /*parent=*/0,
        SE3<S>::Identity());

    // Rear wheel — joint predecessor at rear axle (chassis origin)
    int rear_wheel_idx = builder.add_body(
        Body<S>::with_inertia("rear_wheel", p.wheel_mass,
                               Eigen::Matrix<S,3,1>::Zero(), p.wheel_I),
        RevoluteJoint<S>{Eigen::Matrix<S,3,1>(S(0), S(1), S(0))},
        chassis_idx,
        SE3<S>::Identity());

    // Front wheel — joint predecessor at front axle (wheelbase forward)
    int front_wheel_idx = builder.add_body(
        Body<S>::with_inertia("front_wheel", p.wheel_mass,
                               Eigen::Matrix<S,3,1>::Zero(), p.wheel_I),
        RevoluteJoint<S>{Eigen::Matrix<S,3,1>(S(0), S(1), S(0))},
        chassis_idx,
        SE3<S>{Eigen::Matrix<S,3,3>::Identity(),
               Eigen::Matrix<S,3,1>(p.wheelbase, S(0), S(0))});

    MotorcycleModel<S> mm;
    mm.model            = builder.build();
    mm.rear_wheel_body  = rear_wheel_idx;
    mm.front_wheel_body = front_wheel_idx;
    return mm;
}

// ─── Topology with suspension ─────────────────────────────────────────────────
//
//   world → chassis [FreeJoint]
//           ├── rear_slider  [PrismaticJoint, +Z in chassis]
//           │   └── rear_wheel  [RevoluteJoint, +Y]
//           └── front_slider [PrismaticJoint, +Z in chassis]
//               └── front_wheel [RevoluteJoint, +Y]
//
// q > 0 = fork compressed (wheel moves toward chassis along +Z).
// Spring rest position q0 = 0 (no preload; settles under gravity).
//
template<typename S>
MotorcycleModel<S> build_motorcycle_suspended(const MotorcycleParams<S>& p = {}) {
    ModelBuilder<S> builder;
    Eigen::Matrix<S,3,1> fork_axis{S(0), S(0), S(1)};  // +Z in chassis frame

    // Chassis
    int chassis_idx = builder.add_body(
        Body<S>::with_inertia("chassis", p.chassis_mass, p.chassis_com, p.chassis_I),
        FreeJoint<S>{},
        0, SE3<S>::Identity());

    // ── Rear suspension + wheel ───────────────────────────────────────────────
    // Slider: light unsprung mass representing the rear shock linkage.
    int rear_slider_body = builder.add_body(
        Body<S>::with_inertia("rear_slider", p.fork_mass,
                               Eigen::Matrix<S,3,1>::Zero(),
                               Eigen::Matrix<S,3,3>::Identity() * S(0.01)),
        PrismaticJoint<S>{fork_axis},
        chassis_idx,
        SE3<S>::Identity());   // rear axle at chassis origin
    int rear_fork_joint = rear_slider_body - 1;  // joint index = child_body - 1

    int rear_wheel_idx = builder.add_body(
        Body<S>::with_inertia("rear_wheel", p.wheel_mass,
                               Eigen::Matrix<S,3,1>::Zero(), p.wheel_I),
        RevoluteJoint<S>{Eigen::Matrix<S,3,1>(S(0), S(1), S(0))},
        rear_slider_body,
        SE3<S>::Identity());

    // ── Front suspension + wheel ──────────────────────────────────────────────
    int front_slider_body = builder.add_body(
        Body<S>::with_inertia("front_slider", p.fork_mass,
                               Eigen::Matrix<S,3,1>::Zero(),
                               Eigen::Matrix<S,3,3>::Identity() * S(0.01)),
        PrismaticJoint<S>{fork_axis},
        chassis_idx,
        SE3<S>{Eigen::Matrix<S,3,3>::Identity(),
               Eigen::Matrix<S,3,1>(p.wheelbase, S(0), S(0))});
    int front_fork_joint = front_slider_body - 1;

    int front_wheel_idx = builder.add_body(
        Body<S>::with_inertia("front_wheel", p.wheel_mass,
                               Eigen::Matrix<S,3,1>::Zero(), p.wheel_I),
        RevoluteJoint<S>{Eigen::Matrix<S,3,1>(S(0), S(1), S(0))},
        front_slider_body,
        SE3<S>::Identity());

    MotorcycleModel<S> mm;
    mm.model            = builder.build();
    mm.rear_wheel_body  = rear_wheel_idx;
    mm.front_wheel_body = front_wheel_idx;
    mm.rear_fork_joint  = rear_fork_joint;
    mm.front_fork_joint = front_fork_joint;
    return mm;
}

// ─── Helpers ─────────────────────────────────────────────────────────────────

// Initial state: upright, wheels at ground level (z = wheel_radius),
// all joint angles = 0, velocities = 0.
template<typename S>
State<S> motorcycle_standing_state(const Model<S>& model, S wheel_radius = S(0.3)) {
    State<S> s = model.zero_state();
    // FreeJoint q slice: [qx, qy, qz, qw, px, py, pz]
    // zero_state already sets qw=1; set pz = wheel_radius.
    s.q(6) = wheel_radius;
    return s;
}

// Wire tire contact forces for both wheels into a Simulator.
// Call after building the simulator; the callback runs each step.
template<typename S>
void attach_tire_forces(Simulator<S>&          sim,
                        const MotorcycleModel<S>& mm,
                        const MotorcycleParams<S>& p)
{
    TireContact<S> rear_tire{p.wheel_radius, p.tire_k_n, p.tire_c_n, p.tire_mu};
    TireContact<S> front_tire{p.wheel_radius, p.tire_k_n, p.tire_c_n, p.tire_mu};
    const int ri = mm.rear_wheel_body;
    const int fi = mm.front_wheel_body;

    sim.set_force_callback([=](const Model<S>&, const State<S>&,
                               const KinematicsData<S>& kd)
        -> std::vector<Wrench<S>>
    {
        std::vector<Wrench<S>> f(kd.X_world.size(), Wrench<S>::Zero());
        f[ri] = rear_tire .compute(kd.X_world[ri], kd.v_body[ri]);
        f[fi] = front_tire.compute(kd.X_world[fi], kd.v_body[fi]);
        return f;
    });
}

// Wire suspension springs into a Simulator for the suspended topology.
template<typename S>
void attach_suspension(Simulator<S>&             sim,
                       const MotorcycleModel<S>& mm,
                       const MotorcycleParams<S>& p)
{
    if (mm.rear_fork_joint  >= 0)
        sim.add_spring(mm.rear_fork_joint,
                       SuspensionSpring<S>{p.suspension_k, p.suspension_c, S(0)});
    if (mm.front_fork_joint >= 0)
        sim.add_spring(mm.front_fork_joint,
                       SuspensionSpring<S>{p.suspension_k, p.suspension_c, S(0)});
}

}  // namespace examples
}  // namespace rapide
