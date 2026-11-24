#pragma once
#include "rapide/integration/configuration.hpp"
#include "rapide/model/model.hpp"

namespace rapide::maximal {
inline Eigen::VectorXd velocities(const Model &model, const State &state) {
  Eigen::VectorXd v(model.nv);
  for (size_t i = 0; i < model.bodies.size(); ++i) {
    const int offset = model.bodies[i].velocity_offset;
    if (offset < 0) {
      continue;
    }
    v.segment<3>(offset) = state.bodies[i].omega;
    v.segment<3>(offset + 3) = state.bodies[i].velocity;
  }
  return v;
}

inline void set_velocities(const Model &model, State &state,
                           const Eigen::VectorXd &v) {
  if (v.size() != model.nv || !v.allFinite()) {
    throw std::invalid_argument("Invalid velocity vector");
  }
  for (size_t i = 0; i < model.bodies.size(); ++i) {
    const int offset = model.bodies[i].velocity_offset;
    if (offset < 0) {
      continue;
    }
    state.bodies[i].omega = v.segment<3>(offset);
    state.bodies[i].velocity = v.segment<3>(offset + 3);
  }
}

/** Apply body-axis rotation and world-axis translation increments. Preserve
 *  world angular velocity when correcting a body's orientation. */
inline void retract(const Model &model, State &state,
                    const Eigen::VectorXd &delta) {
  if (delta.size() != model.nv || !delta.allFinite()) {
    throw std::invalid_argument("Invalid pose increment");
  }
  for (size_t i = 0; i < model.bodies.size(); ++i) {
    const int offset = model.bodies[i].velocity_offset;
    if (offset < 0) {
      continue;
    }
    auto &s = state.bodies[i];
    const Eigen::Quaterniond dq = rotation_increment(delta.segment<3>(offset));
    s.pose.orientation =
        integration::apply_body_rotation(s.pose.orientation, dq);
    s.pose.position = integration::advance_linear(
        s.pose.position, delta.segment<3>(offset + 3), 1.);
    s.omega = dq.conjugate() * s.omega;
  }
}

} // namespace rapide::maximal
