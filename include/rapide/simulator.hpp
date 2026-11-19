#pragma once
#include <vector>
#include <functional>
#include "rapide/model/model.hpp"
#include "rapide/algorithm/aba.hpp"
#include "rapide/algorithm/kinematics.hpp"
#include "rapide/integrator/integrator.hpp"
#include "rapide/contact/suspension.hpp"
#include "rapide/math/spatial.hpp"

namespace rapide {

// ForceCallback — called each step; returns external forces per body in WORLD frame.
template<typename S>
using ForceCallback = std::function<
    std::vector<Wrench<S>>(const Model<S>&, const State<S>&, const KinematicsData<S>&)
>;

// TauCallback — called each step; returns full generalised-force vector (nv-dim).
// Its result is ADDED to the constant tau and spring contributions.
template<typename S>
using TauCallback = std::function<Eigen::VectorX<S>(const Model<S>&, const State<S>&)>;

template<typename S>
class Simulator {
public:
    Eigen::Matrix<S,3,1> gravity{S(0), S(0), S(-9.81)};
    S time = S(0);

    explicit Simulator(Model<S> model)
        : model_(std::move(model))
        , state_(model_.zero_state())
        , const_tau_(Eigen::VectorX<S>::Zero(model_.nv))
    {}

    const Model<S>& model()  const { return model_; }
    const State<S>& state()  const { return state_; }
    void set_state(State<S> s)     { state_ = std::move(s); }

    // Constant generalised torques (overwritten each call; not accumulated).
    void set_tau(Eigen::VectorX<S> tau) { const_tau_ = std::move(tau); }

    // Passive spring on a specific joint (accumulated across calls).
    void add_spring(int joint_idx, SuspensionSpring<S> spring) {
        springs_.push_back({joint_idx, std::move(spring)});
    }

    // Active tau callback (replaces any previously set callback).
    void set_tau_callback(TauCallback<S> cb) { tau_cb_ = std::move(cb); }

    // External-force callback (replaces any previously set callback).
    void set_force_callback(ForceCallback<S> cb) { force_cb_ = std::move(cb); }

    // Advance by dt seconds using semi-implicit Euler.
    void step(S dt) {
        KinematicsData<S> kd = forward_kinematics(model_, state_);

        // Assemble external forces.
        std::vector<Wrench<S>> f_ext(model_.num_bodies(), Wrench<S>::Zero());
        if (force_cb_) f_ext = force_cb_(model_, state_, kd);

        // Assemble tau: constant + springs + callback.
        Eigen::VectorX<S> tau = const_tau_;
        apply_springs(springs_, model_, state_, tau);
        if (tau_cb_) tau += tau_cb_(model_, state_);

        Eigen::VectorX<S> ddq = forward_dynamics(model_, state_, tau, f_ext, gravity);
        state_ = euler_step(model_, state_, ddq, dt);
        time += dt;
    }

private:
    Model<S>                    model_;
    State<S>                    state_;
    Eigen::VectorX<S>           const_tau_;
    std::vector<JointSpring<S>> springs_;
    TauCallback<S>              tau_cb_;
    ForceCallback<S>            force_cb_;
};

}  // namespace rapide
