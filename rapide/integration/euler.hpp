#pragma once
#include "rapide/integration/configuration.hpp"
#include <utility>

namespace rapide::integration {

template <typename Configuration, typename Velocity> struct StepResult {
  Configuration configuration;
  Velocity velocity;
};

/** Use the new velocity to advance configuration; acceleration is supplied by
 * dynamics. */
template <typename Q, typename V, typename A, typename S, typename Advance>
inline auto semi_implicit_euler(const Q &configuration, const V &velocity,
                                const A &acceleration, S dt, Advance advance) {
  auto next_velocity = advance_linear(velocity, acceleration, dt);
  auto next_configuration = advance(configuration, next_velocity, dt);
  return StepResult<decltype(next_configuration), decltype(next_velocity)>{
      std::move(next_configuration), std::move(next_velocity)};
}

/** Use the old velocity to advance configuration. */
template <typename Q, typename V, typename A, typename S, typename Advance>
inline auto explicit_euler(const Q &configuration, const V &velocity,
                           const A &acceleration, S dt, Advance advance) {
  auto next_configuration = advance(configuration, velocity, dt);
  auto next_velocity = advance_linear(velocity, acceleration, dt);
  return StepResult<decltype(next_configuration), decltype(next_velocity)>{
      std::move(next_configuration), std::move(next_velocity)};
}

} // namespace rapide::integration
