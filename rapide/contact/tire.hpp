#pragma once
#include "rapide/math/se3.hpp"
#include "rapide/math/spatial.hpp"

namespace rapide {

// TireContact<S> — spring-damper normal force + linear-slip friction
// on a flat ground plane at z = 0.
//
// Returns Wrench in WORLD frame acting on the wheel body.
// Call once per wheel per step; add the result to f_ext[wheel_body_index].
template <typename S> struct TireContact {
  S radius; // unloaded wheel radius
  S k_n;    // normal spring stiffness  [N/m]
  S c_n;    // normal damping           [N·s/m]
  S mu;     // friction coefficient (isotropic, linear slip model)

  // X_wheel_world: SE3 pose of wheel body (centre of wheel) in world frame.
  //   .p = position of wheel centre in world
  //   .R = wheel body orientation
  // v_wheel_world: wheel-body spatial velocity in WHEEL BODY frame.
  Wrench<S> compute(const SE3<S> &X_wheel_world,
                    const Twist<S> &v_wheel_body) const {
    const Eigen::Matrix<S, 3, 1> &p_world = X_wheel_world.p;

    // Penetration depth: how far the wheel centre is below (radius) from
    // ground.
    S penetration = radius - p_world(2);
    if (penetration <= S(0)) {
      return Wrench<S>::Zero(); // no contact
    }

    // Normal velocity at contact patch (world z direction).
    // v_contact_world = v_centre_world (translational part)
    Eigen::Matrix<S, 3, 1> v_centre_world =
        X_wheel_world.R * v_wheel_body.linear;
    S v_normal = v_centre_world(2); // positive = moving up

    // Normal force (spring-damper, clamped to prevent pull).
    S F_n = k_n * penetration - c_n * v_normal;
    if (F_n < S(0)) {
      F_n = S(0);
    }

    // Slip velocity: horizontal component of contact-patch velocity.
    // For a rolling wheel the contact-patch velocity adds the spin
    // contribution:
    //   v_contact = v_centre + ω × r_contact
    // where r_contact is the vector from wheel centre to contact patch (=
    // -radius * z_hat).
    Eigen::Matrix<S, 3, 1> omega_world = X_wheel_world.R * v_wheel_body.angular;
    Eigen::Matrix<S, 3, 1> r_contact(S(0), S(0), -radius);
    Eigen::Matrix<S, 3, 1> v_contact_world =
        v_centre_world + omega_world.cross(r_contact);

    // Horizontal slip velocity (ignore z component).
    Eigen::Matrix<S, 3, 1> v_slip(v_contact_world(0), v_contact_world(1), S(0));

    // Linear Coulomb friction: F_friction = -mu * F_n * normalise(v_slip)
    // (for near-zero slip fall back to zero, avoiding division by zero).
    Eigen::Matrix<S, 3, 1> F_friction = Eigen::Matrix<S, 3, 1>::Zero();
    S slip_mag = v_slip.norm();
    if (slip_mag > S(1e-6)) {
      F_friction = -(mu * F_n / slip_mag) * v_slip;
    }

    // Total force at contact patch in world frame.
    Eigen::Matrix<S, 3, 1> F_total(F_friction(0), F_friction(1), F_n);

    // Torque about wheel centre due to contact force (applied at r_contact).
    Eigen::Matrix<S, 3, 1> tau = r_contact.cross(F_total);

    return {tau, F_total};
  }
};

// Convenience: compute gravitational force on a body in world frame.
template <typename S>
Wrench<S> gravity_force(S mass, const Eigen::Matrix<S, 3, 1> &g = {S(0), S(0),
                                                                   S(-9.81)}) {
  return {Eigen::Matrix<S, 3, 1>::Zero(), mass * g};
}

} // namespace rapide
