#pragma once
#include "rapide/model/forces.hpp"
#include <algorithm>

namespace rapide {

struct TireResponse {
  BodyForce wrench;
  double compression = 0;
  double slip = 0;
  double longitudinal_slip = 0, lateral_slip = 0;
  double longitudinal_force = 0, lateral_force = 0;
  Vec3 contact_point = Vec3::Zero();
  Vec3 slip_velocity = Vec3::Zero();
};

/** Elliptic torus centered at the COM, with local Y as its axle.
 *  Contact uses the lowest support point projected onto a world-space plane
 *  as an approximation to the patch center, without a distributed pressure
 * solve.
 */
struct Tire {
  double radius = .35;
  double rim_radius = .23;
  double width = .14;
  double stiffness = 4000;
  double damping = 100;
  double friction = .8;
  double longitudinal_slip_stiffness = 150;
  double lateral_slip_stiffness = 150;
  double ground_height = 0;
  Vec3 ground_normal = Vec3::UnitZ();

  TireResponse compute(const BodyState &wheel) const {
    if (!std::isfinite(radius) || !std::isfinite(rim_radius) ||
        rim_radius <= 0 || radius <= rim_radius || !std::isfinite(width) ||
        width <= 0 || !std::isfinite(stiffness) || stiffness < 0 ||
        !std::isfinite(damping) || damping < 0 || !std::isfinite(friction) ||
        friction < 0 || !std::isfinite(longitudinal_slip_stiffness) ||
        longitudinal_slip_stiffness < 0 ||
        !std::isfinite(lateral_slip_stiffness) || lateral_slip_stiffness < 0 ||
        !std::isfinite(ground_height) || !ground_normal.allFinite() ||
        std::abs(ground_normal.norm() - 1) > 1e-8) {
      throw std::invalid_argument("Invalid tire parameters");
    }
    const Vec3 axle = wheel.pose.orientation * Vec3::UnitY();
    const Vec3 parallel = ground_normal.dot(axle) * axle;
    const Vec3 radial = ground_normal - parallel;
    const double radial_norm = radial.norm();
    const double major = (radius + rim_radius) / 2;
    const double minor = (radius - rim_radius) / 2;
    const double half_width = width / 2;
    const double support =
        std::sqrt(minor * minor * radial.squaredNorm() +
                  half_width * half_width * parallel.squaredNorm());
    Vec3 support_arm =
        -(minor * minor * radial + half_width * half_width * parallel) /
        support;
    if (radial_norm > 1e-12) {
      support_arm -= major * radial / radial_norm;
    }
    const Vec3 lowest = wheel.pose.position + support_arm;
    const double gap = ground_normal.dot(lowest) - ground_height;
    TireResponse result;
    result.contact_point = lowest - gap * ground_normal;
    result.compression = std::max(0., -gap);
    const Vec3 arm = result.contact_point - wheel.pose.position;
    const Vec3 omega = wheel.pose.orientation * wheel.omega;
    const Vec3 velocity = wheel.velocity + omega.cross(arm);
    result.slip_velocity =
        velocity - ground_normal.dot(velocity) * ground_normal;
    result.slip = result.slip_velocity.norm();
    Vec3 forward = axle.cross(ground_normal);
    if (forward.norm() < 1e-12) {
      forward = wheel.pose.orientation * Vec3::UnitX();
      forward -= ground_normal.dot(forward) * ground_normal;
    }
    forward.normalize();
    const Vec3 lateral = ground_normal.cross(forward);
    result.longitudinal_slip = forward.dot(result.slip_velocity);
    result.lateral_slip = lateral.dot(result.slip_velocity);
    if (result.compression == 0) {
      return result;
    }
    const double normal =
        std::max(0., stiffness * result.compression -
                         damping * ground_normal.dot(velocity));
    Vec3 friction_force =
        -longitudinal_slip_stiffness * result.longitudinal_slip * forward -
        lateral_slip_stiffness * result.lateral_slip * lateral;
    const double magnitude = friction_force.norm();
    if (magnitude > friction * normal) {
      friction_force *= friction * normal / magnitude;
    }
    result.longitudinal_force = forward.dot(friction_force);
    result.lateral_force = lateral.dot(friction_force);
    result.wrench.force = normal * ground_normal + friction_force;
    result.wrench.torque = arm.cross(result.wrench.force);
    return result;
  }
};

using PlanarTire = Tire;

} // namespace rapide
