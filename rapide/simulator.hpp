#pragma once
#include "rapide/aba/backend.hpp"
#include "rapide/maximal/backend.hpp"
#include <functional>
#include <memory>

namespace rapide {

using ForceCallback = std::function<Forces(const Model &, const State &)>;

class Simulator {
public:
  Vec3 gravity{0, 0, -9.81};

  /** Initialize from model.initial_state using the selected backend's policy:
   *  maximal projects joint constraints; ABA requires a compatible tree state.
   */
  explicit Simulator(Model model, Backend backend = Backend::Maximal,
                     SolverOptions options = {})
      : Simulator(std::move(model), make_backend(backend), options) {}

  Simulator(Model model, SolverOptions options)
      : Simulator(std::move(model), Backend::Maximal, options) {}

  Simulator(Model model, std::unique_ptr<SimulationBackend> backend,
            SolverOptions options = {})
      : model_(std::move(model)), backend_(std::move(backend)) {
    if (!backend_) {
      throw std::invalid_argument("Simulator needs a backend");
    }
    model_.validate();
    backend_->setup(model_, options);
  }

  const Model &model() const { return model_; }
  const State &state() const { return backend_->state(); }
  const StepReport &report() const { return backend_->report(); }
  double time() const { return time_; }
  /** Replace body state using the backend's initialization policy. Simulation
   *  time is unchanged; failed assignments preserve the current state. */
  void set_state(State state) { backend_->set_state(std::move(state)); }
  void set_force_callback(ForceCallback callback) {
    callback_ = std::move(callback);
  }

  void set_joint_effort(JointId joint, double effort) {
    joint_position(model_, state(), joint);
    if (!std::isfinite(effort)) {
      throw std::invalid_argument("Nonfinite effort");
    }
    if (efforts_.empty()) {
      efforts_.resize(model_.joints.size());
    }
    efforts_[joint] = effort;
  }

  void add_joint_spring(JointId joint, double stiffness, double damping,
                        double rest) {
    auto check = zero_forces(model_);
    rapide::add_joint_spring(model_, state(), check, joint, stiffness, damping,
                             rest);
    springs_.push_back({joint, stiffness, damping, rest});
  }

  StepReport step(double dt) {
    auto forces = callback_ ? callback_(model_, state()) : zero_forces(model_);
    if (forces.size() != model_.bodies.size()) {
      throw std::invalid_argument("Force count does not match model");
    }
    for (size_t i = 0; i < efforts_.size(); ++i) {
      if (efforts_[i] != 0) {
        add_joint_effort(model_, state(), forces, static_cast<int>(i),
                         efforts_[i]);
      }
    }
    for (const auto &spring : springs_) {
      rapide::add_joint_spring(model_, state(), forces, spring.joint,
                               spring.stiffness, spring.damping, spring.rest);
    }
    return step(dt, forces);
  }

  StepReport step(double dt, const Forces &forces) {
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
    const auto report = backend_->advance(dt, forces, gravity);
    time_ += dt;
    return report;
  }

private:
  static std::unique_ptr<SimulationBackend> make_backend(Backend backend) {
    switch (backend) {
    case Backend::ABA:
      return std::make_unique<aba::ABABackend>();
    case Backend::Maximal:
      return std::make_unique<maximal::MaximalBackend>();
    }
    throw std::invalid_argument("Unknown backend");
  }
  struct JointSpring {
    JointId joint;
    double stiffness, damping, rest;
  };
  std::vector<JointSpring> springs_;
  std::vector<double> efforts_;
  Model model_;
  std::unique_ptr<SimulationBackend> backend_;
  ForceCallback callback_;
  double time_ = 0;
};

} // namespace rapide
