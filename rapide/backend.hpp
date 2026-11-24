#pragma once
#include "rapide/model/forces.hpp"

namespace rapide {

enum class Backend { ABA, Maximal };

struct SolverOptions {
  double position_tolerance = 1e-9;
  double orientation_tolerance = 1e-9;
  double velocity_tolerance = 1e-8;
  double rank_threshold = 1e-10;
  double max_translation_correction = .1;
  double max_rotation_correction = .25;
  int max_projection_iterations = 30;

  void validate() const {
    for (double x : {position_tolerance, orientation_tolerance,
                     velocity_tolerance, rank_threshold,
                     max_translation_correction, max_rotation_correction}) {
      if (!std::isfinite(x) || x <= 0) {
        throw std::invalid_argument(
            "Solver tolerances must be positive and finite");
      }
    }
    if (rank_threshold >= 1 || max_projection_iterations < 1) {
      throw std::invalid_argument("Invalid rank threshold or iteration limit");
    }
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

/** Backends compile a shared model and expose shared COM body states. */
class SimulationBackend {
public:
  virtual ~SimulationBackend() = default;
  virtual void setup(const Model &model, const SolverOptions &options) = 0;
  virtual const State &state() const = 0;
  virtual const StepReport &report() const = 0;
  virtual void set_state(State state) = 0;
  virtual StepReport advance(double dt, const Forces &forces,
                             const Vec3 &gravity) = 0;
};

} // namespace rapide
