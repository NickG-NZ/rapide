#pragma once
#include <vector>
#include "rapide/model/model.hpp"

namespace rapide {

// Linear spring-damper for a single passive joint DOF.
// Attach to any 1-DOF joint (typically PrismaticJoint for fork suspension,
// or RevoluteJoint for a swing-arm pivot).
template<typename S>
struct SuspensionSpring {
    S k;   // stiffness  [N/m for prismatic, N·m/rad for revolute]
    S c;   // damping    [N·s/m or N·m·s/rad]
    S q0;  // rest (zero-force) position

    // Generalised force contribution for this joint's DOF.
    S compute_tau(S q, S v) const {
        return -k * (q - q0) - c * v;
    }
};

// A spring attached to a specific joint by index.
template<typename S>
struct JointSpring {
    int                 joint_idx;
    SuspensionSpring<S> spring;
};

// Compute all spring tau contributions and add them to tau in-place.
template<typename S>
void apply_springs(const std::vector<JointSpring<S>>& springs,
                   const Model<S>&                    model,
                   const State<S>&                    state,
                   Eigen::VectorX<S>&                 tau)
{
    for (const auto& js : springs) {
        const int qi = model.q_idx[js.joint_idx];
        const int vi = model.v_idx[js.joint_idx];
        tau(vi) += js.spring.compute_tau(state.q(qi), state.v(vi));
    }
}

}  // namespace rapide
