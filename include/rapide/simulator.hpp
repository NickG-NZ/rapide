#pragma once
#include <vector>
#include <functional>
#include "rapide/model/model.hpp"
#include "rapide/algorithm/aba.hpp"
#include "rapide/algorithm/kinematics.hpp"
#include "rapide/integrator/integrator.hpp"
#include "rapide/math/spatial.hpp"

namespace rapide {

// ForceCallback: called each step with kinematics data; returns f_ext per body.
// Use this to wire in tire forces, controller forces, or any body-frame loads.
// Forces should be expressed in WORLD frame.
template<typename S>
using ForceCallback = std::function<
    std::vector<Force6<S>>(const Model<S>&, const State<S>&, const KinematicsData<S>&)
>;

template<typename S>
class Simulator {
public:
    Eigen::Matrix<S,3,1> gravity{S(0), S(0), S(-9.81)};
    S time = S(0);

    explicit Simulator(Model<S> model)
        : model_(std::move(model))
        , state_(model_.zero_state())
    {}

    const Model<S>& model()  const { return model_; }
    const State<S>& state()  const { return state_; }
    void set_state(State<S> s)     { state_ = std::move(s); }

    // Set the control input (generalised torques, nv-dimensional).
    void set_tau(Eigen::VectorX<S> tau) { tau_ = std::move(tau); }

    // Register a callback that computes external forces each step.
    void set_force_callback(ForceCallback<S> cb) { force_cb_ = std::move(cb); }

    // Advance the simulation by dt seconds (semi-implicit Euler).
    void step(S dt) {
        // Forward kinematics (needed for force callbacks and contact).
        KinematicsData<S> kd = forward_kinematics(model_, state_);

        // Assemble external forces (zeros by default, filled by callback).
        std::vector<Force6<S>> f_ext(model_.num_bodies(), Force6<S>::Zero());
        if (force_cb_) {
            f_ext = force_cb_(model_, state_, kd);
        }

        // Generalised accelerations via ABA.
        Eigen::VectorX<S> tau = tau_.size() == model_.nv
            ? tau_
            : Eigen::VectorX<S>::Zero(model_.nv);

        Eigen::VectorX<S> ddq = forward_dynamics(model_, state_, tau, f_ext, gravity);

        // Integrate state.
        state_ = euler_step(model_, state_, ddq, dt);
        time += dt;
    }

private:
    Model<S>         model_;
    State<S>         state_;
    Eigen::VectorX<S> tau_;
    ForceCallback<S> force_cb_;
};

}  // namespace rapide
