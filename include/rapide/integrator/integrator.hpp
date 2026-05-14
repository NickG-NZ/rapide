#pragma once
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "rapide/model/model.hpp"
#include "rapide/model/joint.hpp"

namespace rapide {

// Integrate joint configurations forward by dt given current velocities.
// Handles per-joint kinematics (quaternion integration for FreeJoint).
template<typename S>
Eigen::VectorX<S> integrate_configuration(
    const Model<S>&          model,
    const Eigen::VectorX<S>& q,
    const Eigen::VectorX<S>& v,
    S dt)
{
    Eigen::VectorX<S> q_new = q;
    for (int i = 0; i < model.num_joints(); ++i) {
        std::visit([&](const auto& jt) {
            using JT = std::decay_t<decltype(jt)>;
            if constexpr (std::is_same_v<JT, FreeJoint<S>>) {
                Eigen::Matrix<S,7,1> qi = model.q_i(q, i);
                Eigen::Matrix<S,6,1> vi = model.v_i(v, i);
                q_new.segment(model.q_idx[i], 7) = jt.integrate_q(qi, vi, dt);
            } else {
                // RevoluteJoint: trivial q += v * dt
                q_new(model.q_idx[i]) += v(model.v_idx[i]) * dt;
            }
        }, model.joints[i]);
    }
    return q_new;
}

// Semi-implicit (symplectic) Euler integrator.
//
// Update order: v_new = v + ddq * dt,  q_new = integrate(q, v_new, dt).
// More stable than explicit Euler for conservative mechanical systems.
template<typename S>
State<S> euler_step(const Model<S>& model,
                    const State<S>& state,
                    const Eigen::VectorX<S>& ddq,
                    S dt)
{
    Eigen::VectorX<S> v_new = state.v + ddq * dt;
    Eigen::VectorX<S> q_new = integrate_configuration(model, state.q, v_new, dt);
    return {q_new, v_new};
}

// Explicit Euler (for comparison / debugging).
template<typename S>
State<S> euler_step_explicit(const Model<S>& model,
                              const State<S>& state,
                              const Eigen::VectorX<S>& ddq,
                              S dt)
{
    Eigen::VectorX<S> q_new = integrate_configuration(model, state.q, state.v, dt);
    Eigen::VectorX<S> v_new = state.v + ddq * dt;
    return {q_new, v_new};
}

}  // namespace rapide
