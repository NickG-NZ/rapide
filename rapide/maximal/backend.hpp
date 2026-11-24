#pragma once
#include "rapide/backend.hpp"
#include "rapide/maximal/constraints.hpp"
#include <Eigen/QR>
#include <functional>
#include <utility>

namespace rapide::maximal {

/** Dense, rank-aware bilateral impulse solver with mass-weighted position
 *  projection. Steps and state assembly commit only after convergence. */
class MaximalBackend final : public SimulationBackend {
public:
  void setup(const Model &model, const SolverOptions &options) override {
    model_ = model;
    state_ = model.initial_state;
    options_ = options;
    options_.validate();
    model_.validate_state(state_);
    report_ = assemble(state_);
  }

  const State &state() const override { return state_; }
  const StepReport &report() const override { return report_; }

  /** Locally project supplied poses and velocities onto the bilateral joints.
   *  Throws if assembly fails; the previously accepted state is preserved.
   *  This enforces joint compatibility, not spring/contact equilibrium. */
  void set_state(State state) override {
    model_.validate_state(state);
    auto report = assemble(state);
    state_ = std::move(state);
    report_ = std::move(report);
  }

  /** Predict velocities from world forces/COM torques and gyroscopic bias,
   *  then apply joint impulses to enforce J*v = 0. Integrate poses using those
   *  velocities and project away position and velocity drift before committing.
   *  Loads stay constant over dt. Projection impulses are reported separately
   *  and are not physical reaction-force measurements. */
  StepReport advance(double dt, const Forces &forces,
                     const Vec3 &gravity) override {
    if (!std::isfinite(dt) || dt <= 0 || !gravity.allFinite()) {
      throw std::invalid_argument("Invalid timestep or gravity");
    }
    if (forces.size() != model_.bodies.size()) {
      throw std::invalid_argument("Force count does not match model");
    }
    for (const auto &force : forces) {
      if (!force.force.allFinite() || !force.torque.allFinite()) {
        throw std::invalid_argument("Nonfinite force");
      }
    }
    State next = state_;
    Eigen::VectorXd v = velocities(model_, next);
    for (size_t i = 0; i < model_.bodies.size(); ++i) {
      const auto &body = model_.bodies[i];
      if (!body.dynamic()) {
        continue;
      }
      const auto &s = next.bodies[i];
      const Vec3 angular_acceleration =
          body.inverse_inertia *
          (s.pose.orientation.conjugate() * forces[i].torque -
           s.omega.cross(body.inertia * s.omega));
      const Vec3 linear_acceleration = gravity + forces[i].force / body.mass;
      v.segment<3>(body.velocity_offset) =
          integration::advance_linear(s.omega, angular_acceleration, dt);
      v.segment<3>(body.velocity_offset + 3) =
          integration::advance_linear(s.velocity, linear_acceleration, dt);
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
      if (!body.dynamic()) {
        continue;
      }
      const int o = body.velocity_offset;
      out.middleRows(o, 3) = body.inverse_inertia * rhs.middleRows(o, 3);
      out.middleRows(o + 3, 3) = rhs.middleRows(o + 3, 3) / body.mass;
    }
    return out;
  }

  /** Solve (J*M^-1*J^T)*p = rhs and return delta = M^-1*J^T*p.
   *  For consistent rows, delta minimizes 0.5*delta^T*M*delta subject to
   *  J*delta = rhs. M uses body-frame inertia and world-frame COM mass blocks.
   *  Diagonal scaling and rank-revealing decomposition handle redundant rows;
   *  incompatible rows produce a least-squares result checked by assembly. */
  Solution solve(const Eigen::MatrixXd &j, const Eigen::VectorXd &rhs) const {
    if (j.rows() == 0) {
      return {Eigen::VectorXd::Zero(model_.nv), Eigen::VectorXd{}, 0};
    }
    const Eigen::MatrixXd response = inverse_mass_times(j.transpose());
    const Eigen::MatrixXd w = j * response;
    Eigen::VectorXd scale(w.rows());
    for (int i = 0; i < w.rows(); ++i) {
      scale(i) = w(i, i) > 0 ? 1 / std::sqrt(w(i, i)) : 1;
    }
    const Eigen::MatrixXd normalized =
        scale.asDiagonal() * w * scale.asDiagonal();
    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> decomposition;
    decomposition.setThreshold(options_.rank_threshold);
    decomposition.compute(normalized);
    const Eigen::VectorXd p =
        scale.asDiagonal() * decomposition.solve(scale.asDiagonal() * rhs);
    const Eigen::VectorXd delta = response * p;
    if (!delta.allFinite() || !p.allFinite()) {
      throw std::runtime_error("Nonfinite constraint solve");
    }
    return {delta, p, static_cast<int>(decomposition.rank())};
  }

  double error(const Constraints &c) const {
    return std::max(c.position_error() / options_.position_tolerance,
                    c.orientation_error() / options_.orientation_tolerance);
  }

  /** Assemble nonlinear joint constraints with repeated mass-weighted solves
   *  J*delta = -residual, quaternion retraction, and backtracking for descent.
   *  Once poses converge, project velocities onto J*v = 0 at the final poses.
   *  This is a local feasibility solve, not a globally nearest-pose optimizer;
   *  inconsistent geometry or a poor initial guess can prevent convergence. */
  StepReport assemble(State &state) const {
    StepReport report;
    auto constraints = assemble_constraints(model_, state);
    for (; error(constraints) > 1; ++report.projection_iterations) {
      if (report.projection_iterations >= options_.max_projection_iterations) {
        throw std::runtime_error("Position projection did not converge; "
                                 "inspect linkage geometry and singularities");
      }
      const auto correction =
          solve(constraints.jacobian, -constraints.residual);
      double alpha = 1;
      for (const auto &body : model_.bodies) {
        if (!body.dynamic()) {
          continue;
        }
        const int o = body.velocity_offset;
        const double rotation = correction.delta.segment<3>(o).norm();
        const double translation = correction.delta.segment<3>(o + 3).norm();
        if (rotation > 0) {
          alpha = std::min(alpha, options_.max_rotation_correction / rotation);
        }
        if (translation > 0) {
          alpha = std::min(alpha,
                           options_.max_translation_correction / translation);
        }
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
      if (!accepted) {
        throw std::runtime_error(
            "Position projection stalled; constraints may be inconsistent");
      }
    }
    Eigen::VectorXd v = velocities(model_, state);
    const auto correction =
        solve(constraints.jacobian, -constraints.jacobian * v);
    v += correction.delta;
    const Eigen::VectorXd residual = constraints.jacobian * v;
    report.velocity_error =
        residual.size() ? residual.lpNorm<Eigen::Infinity>() : 0;
    if (report.velocity_error > options_.velocity_tolerance) {
      throw std::runtime_error("Velocity projection did not converge");
    }
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
};

} // namespace rapide::maximal
