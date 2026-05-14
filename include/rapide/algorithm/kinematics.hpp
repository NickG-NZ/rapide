#pragma once
#include <vector>
#include "rapide/model/model.hpp"
#include "rapide/math/se3.hpp"

namespace rapide {

// Forward kinematics: compute SE3 pose of every body in world frame.
//
// Returns X_world[i] for i in [0, num_bodies):
//   X_world[0] = identity (world)
//   X_world[i] = pose of bodies[i] in world frame for i >= 1
//
// Also fills body-frame velocities v_body[i] (Motion6 in body frame).
template<typename S>
struct KinematicsData {
    std::vector<SE3<S>>     X_world;   // body pose in world frame
    std::vector<Motion6<S>> v_body;    // spatial velocity in body frame
};

template<typename S>
KinematicsData<S> forward_kinematics(const Model<S>& model, const State<S>& state) {
    const int nb = model.num_bodies();
    KinematicsData<S> kd;
    kd.X_world.resize(nb, SE3<S>::Identity());
    kd.v_body.resize(nb, Motion6<S>::Zero());

    for (int i = 0; i < model.num_joints(); ++i) {
        const int child_body = i + 1;
        const int par        = model.parent[i];

        // Compute X_J: transform across joint from predecessor to successor frame.
        SE3<S> X_J = std::visit([&](const auto& jt) -> SE3<S> {
            using JT = std::decay_t<decltype(jt)>;
            if constexpr (std::is_same_v<JT, FreeJoint<S>>) {
                Eigen::Matrix<S,7,1> qi = model.q_i(state.q, i);
                return jt.X_J(qi);
            } else {
                return jt.X_J(state.q(model.q_idx[i]));
            }
        }, model.joints[i]);

        // Body pose in world = parent pose * fixed tree transform * joint transform.
        kd.X_world[child_body] = kd.X_world[par] * model.X_tree[i] * X_J;

        // Joint velocity contribution in body (successor) frame.
        Motion6<S> vJ = std::visit([&](const auto& jt) -> Motion6<S> {
            using JT = std::decay_t<decltype(jt)>;
            if constexpr (std::is_same_v<JT, FreeJoint<S>>) {
                Eigen::Matrix<S,6,1> vi = model.v_i(state.v, i);
                return jt.vJ(vi);
            } else {
                return jt.vJ(state.v(model.v_idx[i]));
            }
        }, model.joints[i]);

        // Body velocity = transform parent velocity into child frame + joint velocity.
        // v_child = X_J⁻¹ * X_tree⁻¹ * v_parent_in_world ... expressed in body frame.
        // More directly: v_child = (X_tree * X_J).inv_act_on_motion(v_world_parent) + vJ
        SE3<S> X_up = model.X_tree[i] * X_J;  // from parent body to child body
        kd.v_body[child_body] = X_up.inv_act_on_motion(kd.v_body[par]) + vJ;
    }

    return kd;
}

}  // namespace rapide
