#pragma once
#include "rapide/model/model.hpp"
#include "rapide/model/joint.hpp"
#include "rapide/model/body.hpp"
#include "rapide/math/se3.hpp"

namespace rapide {
namespace examples {

// build_arm<S>() — construct an n-DOF serial revolute-joint arm.
//
// Attaches to parent_body_idx in an existing model (typically the chassis).
// Joint axes alternate: Z, Y, Z, ... (a common DH-style convention).
// Link lengths are uniform (link_length).
//
// Returns the modified ModelBuilder with the arm attached.
template<typename S>
void attach_arm(ModelBuilder<S>& builder,
                int parent_body_idx,
                int n_dof,
                S link_length = S(0.4),
                S link_mass   = S(2.0),
                const SE3<S>& X_mount = SE3<S>::Identity())
{
    // Joint axes: Z, Y, Z, ... (shoulder yaw, shoulder pitch, elbow pitch, ...)
    const Eigen::Matrix<S,3,1> axes[3] = {
        {S(0), S(0), S(1)},  // Z
        {S(0), S(1), S(0)},  // Y
        {S(0), S(0), S(1)},  // Z
    };

    Eigen::Matrix<S,3,3> I_link;
    I_link << S(0.01), S(0), S(0),
              S(0), S(0.01), S(0),
              S(0), S(0), S(0.001);

    int parent = parent_body_idx;
    SE3<S> X_parent = X_mount;

    for (int j = 0; j < n_dof; ++j) {
        Body<S> link = Body<S>::with_inertia(
            "arm_link_" + std::to_string(j),
            link_mass,
            Eigen::Matrix<S,3,1>(link_length / S(2), S(0), S(0)),  // CoM at mid-link
            I_link);

        Eigen::Matrix<S,3,1> axis = axes[j % 3];
        parent = builder.add_body(
            std::move(link),
            RevoluteJoint<S>{axis},
            parent,
            X_parent);

        // Next joint is at the distal end of this link (+X in link frame).
        X_parent = SE3<S>{
            Eigen::Matrix<S,3,3>::Identity(),
            Eigen::Matrix<S,3,1>(link_length, S(0), S(0))
        };
    }

    // End-effector mass (simulates human hand / payload).
    Body<S> ee = Body<S>::point_mass(
        "arm_ee",
        S(5.0),
        Eigen::Matrix<S,3,1>::Zero());
    builder.add_body(
        std::move(ee),
        RevoluteJoint<S>{{S(0), S(1), S(0)}},  // wrist (placeholder — locked by τ=0)
        parent,
        X_parent);
}

}  // namespace examples
}  // namespace rapide
