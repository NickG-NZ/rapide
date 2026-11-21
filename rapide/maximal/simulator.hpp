#pragma once
#include "rapide/maximal/constraints.hpp"
#include "rapide/maximal/forces.hpp"
#include <Eigen/QR>
#include <functional>
#include <utility>

namespace rapide::maximal {

struct SolverOptions {
  double position_tolerance = 1e-9;
  double orientation_tolerance = 1e-9;
  double velocity_tolerance = 1e-8;
  double rank_threshold = 1e-10;
  double max_translation_correction = .1;
  double max_rotation_correction = .25;
  int max_projection_iterations = 30;

  void validate() const {
    for (double x :
         {position_tolerance, orientation_tolerance, velocity_tolerance,
          rank_threshold, max_translation_correction, max_rotation_correction})
      if (!std::isfinite(x) || x <= 0)
        throw std::invalid_argument(
            "Solver tolerances must be positive and finite");
    if (rank_threshold >= 1 || max_projection_iterations < 1)
      throw std::invalid_argument("Invalid rank threshold or iteration limit");
  }
};

struct StepReport {
  int constraint_rows = 0;
  int rank = 0;
  int projection_iterations = 0;
  double position_error = 0;
  double orientation_error = 0;
  double velocity_error = 0;
  Eigen::VectorXd constraint_impulses;
  Eigen::VectorXd projection_impulses;
};

using ForceCallback = std::function<Forces(const Model &, const State &)>;

/** Dense, rank-aware bilateral impulse solver with mass-weighted position
 *  projection. Steps and state assembly commit only after convergence. */
class Simulator {
public:
  Vec3 gravity{0, 0, -9.81};

  explicit Simulator(Model model, SolverOptions options = {})
      : model_(std::move(model)), state_(model_.initial_state),
        options_(options) {
    options_.validate();
    model_.validate_state(state_);
    report_ = assemble(state_);
  }

  const Model &model() const { return model_; }
  const State &state() const { return state_; }
  const StepReport &report() const { return report_; }
  double time() const { return time_; }
  void set_force_callback(ForceCallback callback) {
    callback_ = std::move(callback);
  }

  /** Project a supplied configuration and velocity onto the joint constraints.
   */
  void set_state(State state) {
    model_.validate_state(state);
    auto report = assemble(state);
    state_ = std::move(state);
    report_ = std::move(report);
  }

  StepReport step(double dt) {
    return step(dt,
                callback_ ? callback_(model_, state_) : zero_forces(model_));
  }

  /** Advance with world forces/COM torques held constant over dt. Constraint
   *  impulses precede pose integration; projection impulses are reported
   *  separately and are not physical reaction-force measurements. */
  StepReport step(double dt, const Forces &forces) {
    if (!std::isfinite(dt) || dt <= 0 || !gravity.allFinite())
      throw std::invalid_argument("Invalid timestep or gravity");
    if (forces.size() != model_.bodies.size())
      throw std::invalid_argument("Force count does not match model");
    for (const auto &force : forces)
      if (!force.force.allFinite() || !force.torque.allFinite())
        throw std::invalid_argument("Nonfinite force");
    State next = state_;
    Eigen::VectorXd v = velocities(model_, next);
    for (size_t i = 0; i < model_.bodies.size(); ++i) {
      const auto &body = model_.bodies[i];
      if (!body.dynamic())
        continue;
      const auto &s = next.bodies[i];
      v.segment<3>(body.velocity_offset) +=
          dt * body.inverse_inertia *
          (s.pose.orientation.conjugate() * forces[i].torque -
           s.omega.cross(body.inertia * s.omega));
      v.segment<3>(body.velocity_offset + 3) +=
          dt * (gravity + forces[i].force / body.mass);
    }
    auto constraints = assemble_constraints(model_, next);
    const auto impulse = solve(constraints.jacobian, -constraints.jacobian * v);
    v += impulse.delta;
    set_velocities(model_, next, v);
    retract(model_, next, dt * v);
    StepReport report = assemble(next);
    report.constraint_impulses = impulse.impulses;
    model_.validate_state(next);
    state_ = std::move(next);
    report_ = std::move(report);
    time_ += dt;
    return report_;
  }

private:
  struct Solution {
    Eigen::VectorXd delta, impulses;
    int rank;
  };

  Eigen::MatrixXd inverse_mass_times(const Eigen::MatrixXd &rhs) const {
    Eigen::MatrixXd out(rhs.rows(), rhs.cols());
    for (const auto &body : model_.bodies) {
      if (!body.dynamic())
        continue;
      const int o = body.velocity_offset;
      out.middleRows(o, 3) = body.inverse_inertia * rhs.middleRows(o, 3);
      out.middleRows(o + 3, 3) = rhs.middleRows(o + 3, 3) / body.mass;
    }
    return out;
  }

  Solution solve(const Eigen::MatrixXd &j, const Eigen::VectorXd &rhs) const {
    if (j.rows() == 0)
      return {Eigen::VectorXd::Zero(model_.nv), Eigen::VectorXd{}, 0};
    const Eigen::MatrixXd response = inverse_mass_times(j.transpose());
    const Eigen::MatrixXd w = j * response;
    Eigen::VectorXd scale(w.rows());
    for (int i = 0; i < w.rows(); ++i)
      scale(i) = w(i, i) > 0 ? 1 / std::sqrt(w(i, i)) : 1;
    const Eigen::MatrixXd normalized =
        scale.asDiagonal() * w * scale.asDiagonal();
    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> decomposition;
    decomposition.setThreshold(options_.rank_threshold);
    decomposition.compute(normalized);
    const Eigen::VectorXd p =
        scale.asDiagonal() * decomposition.solve(scale.asDiagonal() * rhs);
    const Eigen::VectorXd delta = response * p;
    if (!delta.allFinite() || !p.allFinite())
      throw std::runtime_error("Nonfinite constraint solve");
    return {delta, p, static_cast<int>(decomposition.rank())};
  }

  double error(const Constraints &c) const {
    return std::max(c.position_error() / options_.position_tolerance,
                    c.orientation_error() / options_.orientation_tolerance);
  }

  StepReport assemble(State &state) const {
    StepReport report;
    auto constraints = assemble_constraints(model_, state);
    for (; error(constraints) > 1; ++report.projection_iterations) {
      if (report.projection_iterations >= options_.max_projection_iterations)
        throw std::runtime_error("Position projection did not converge; "
                                 "inspect linkage geometry and singularities");
      const auto correction =
          solve(constraints.jacobian, -constraints.residual);
      double alpha = 1;
      for (const auto &body : model_.bodies) {
        if (!body.dynamic())
          continue;
        const int o = body.velocity_offset;
        const double rotation = correction.delta.segment<3>(o).norm();
        const double translation = correction.delta.segment<3>(o + 3).norm();
        if (rotation > 0)
          alpha = std::min(alpha, options_.max_rotation_correction / rotation);
        if (translation > 0)
          alpha = std::min(alpha,
                           options_.max_translation_correction / translation);
      }
      bool accepted = false;
      for (int attempt = 0; attempt < 16; ++attempt, alpha *= .5) {
        State candidate = state;
        retract(model_, candidate, alpha * correction.delta);
        auto trial = assemble_constraints(model_, candidate);
        if (error(trial) < error(constraints)) {
          state = std::move(candidate);
          constraints = std::move(trial);
          accepted = true;
          break;
        }
      }
      if (!accepted)
        throw std::runtime_error(
            "Position projection stalled; constraints may be inconsistent");
    }
    Eigen::VectorXd v = velocities(model_, state);
    const auto correction =
        solve(constraints.jacobian, -constraints.jacobian * v);
    v += correction.delta;
    const Eigen::VectorXd residual = constraints.jacobian * v;
    report.velocity_error =
        residual.size() ? residual.lpNorm<Eigen::Infinity>() : 0;
    if (report.velocity_error > options_.velocity_tolerance)
      throw std::runtime_error("Velocity projection did not converge");
    set_velocities(model_, state, v);
    report.rank = correction.rank;
    report.constraint_rows = static_cast<int>(constraints.residual.size());
    report.position_error = constraints.position_error();
    report.orientation_error = constraints.orientation_error();
    report.projection_impulses = correction.impulses;
    return report;
  }

  Model model_;
  State state_;
  SolverOptions options_;
  StepReport report_;
  ForceCallback callback_;
  double time_ = 0;
};

} // namespace rapide::maximal
