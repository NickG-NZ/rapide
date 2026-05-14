#pragma once
#include <vector>
#include <variant>
#include <Eigen/Core>
#include <Eigen/LU>
#include "rapide/model/model.hpp"
#include "rapide/math/spatial.hpp"
#include "rapide/math/se3.hpp"

namespace rapide {

// forward_dynamics — Featherstone Articulated Body Algorithm (ABA).
//
// Returns ddq: generalised accelerations (nv-dimensional).
//
// f_ext: external spatial forces per body in WORLD frame (indexed 0..nb-1).
//        Use Force6::Zero() for bodies with no external loads.
//        Gravity is handled internally via the base-acceleration trick;
//        do NOT include gravity in f_ext.
//
// gravity: world-frame gravitational acceleration, e.g. {0, 0, -9.81}.
template<typename S>
Eigen::VectorX<S> forward_dynamics(
    const Model<S>&              model,
    const State<S>&              state,
    const Eigen::VectorX<S>&     tau,
    const std::vector<Force6<S>>& f_ext,
    const Eigen::Matrix<S,3,1>&  gravity = Eigen::Matrix<S,3,1>(S(0), S(0), S(-9.81)))
{
    const int nb = model.num_bodies();
    const int nj = model.num_joints();

    // ── Workspace ────────────────────────────────────────────────────────────
    std::vector<SE3<S>>                X_up(nj);                   // parent→child (see below)
    std::vector<SE3<S>>                X_world(nb, SE3<S>::Identity());
    std::vector<Motion6<S>>            v(nb, Motion6<S>::Zero());
    std::vector<Motion6<S>>            c(nj, Motion6<S>::Zero());
    std::vector<ArticulatedInertia<S>> IA(nb, ArticulatedInertia<S>::Zero());
    std::vector<Force6<S>>             pA(nb, Force6<S>::Zero());
    std::vector<Motion6<S>>            a(nb, Motion6<S>::Zero());

    // Per-joint quantities retained from pass 2 for use in pass 3.
    // Stored as dynamic Eigen types to accommodate 1-DOF and 6-DOF joints.
    std::vector<Eigen::MatrixX<S>>  D_inv(nj);
    std::vector<Eigen::VectorX<S>>  u_vec(nj);
    std::vector<Eigen::MatrixX<S>>  U_mat(nj);  // 6×nv_i

    // Gravity trick: set world-body acceleration to -g so gravity propagates
    // automatically through the recursion (RBDA Section 9.3).
    a[0] = {Eigen::Matrix<S,3,1>::Zero(), -gravity};

    // ── Pass 1: kinematics ───────────────────────────────────────────────────
    // Traverse joints in parent-first order (guaranteed by model construction).
    for (int i = 0; i < nj; ++i) {
        const int bi = i + 1;           // child body index
        const int pi = model.parent[i]; // parent body index

        // Joint transform X_J and joint velocity vJ, dispatched by joint type.
        SE3<S>     X_J;
        Motion6<S> vJ;
        std::visit([&](const auto& jt) {
            using JT = std::decay_t<decltype(jt)>;
            if constexpr (std::is_same_v<JT, FreeJoint<S>>) {
                Eigen::Matrix<S,7,1> qi = model.q_i(state.q, i);
                Eigen::Matrix<S,6,1> vi = model.v_i(state.v, i);
                X_J = jt.X_J(qi);
                vJ  = jt.vJ(vi);
            } else {
                X_J = jt.X_J(state.q(model.q_idx[i]));
                vJ  = jt.vJ(state.v(model.v_idx[i]));
            }
        }, model.joints[i]);

        // X_up[i]: combined fixed+joint transform, child-relative-to-parent.
        // act_on_motion: child→parent; inv_act_on_motion: parent→child.
        X_up[i] = model.X_tree[i] * X_J;

        // Accumulate body pose in world frame.
        X_world[bi] = X_world[pi] * X_up[i];

        // Body velocity in body frame: transform parent velocity into child frame + joint vel.
        v[bi] = X_up[i].inv_act_on_motion(v[pi]) + vJ;

        // Coriolis bias acceleration: c_i = v_i × vJ_i (motion cross product).
        c[i] = v[bi].cross(vJ);

        // Initialise articulated inertia from rigid-body inertia.
        IA[bi] = ArticulatedInertia<S>::from_rigid(model.bodies[bi].inertia);

        // Bias force: p_A = v × I·v  (centrifugal/Coriolis)  minus external force.
        // f_ext is in world frame; convert to body frame via inv_act_on_force.
        Force6<S> fext_body = X_world[bi].inv_act_on_force(f_ext[bi]);
        pA[bi] = cross_motion_force(v[bi], model.bodies[bi].inertia.apply(v[bi])) - fext_body;
    }

    // ── Pass 2: propagate articulated inertia and bias from tips to root ─────
    for (int i = nj - 1; i >= 0; --i) {
        const int bi = i + 1;
        const int pi = model.parent[i];


        // Build U = I_A · S  and  D = Sᵀ · U, dispatched by joint type.
        // U has shape (6 × nv_i), D has shape (nv_i × nv_i).
        std::visit([&](const auto& jt) {
            using JT = std::decay_t<decltype(jt)>;

            if constexpr (std::is_same_v<JT, RevoluteJoint<S>>) {
                // S = [axis; 0] — single column.
                Motion6<S> S_col{jt.axis, Eigen::Matrix<S,3,1>::Zero()};
                Force6<S> U_force = IA[bi].apply(S_col);
                // Store U as 6×1
                Eigen::Matrix<S,6,1> U_vec6;
                U_vec6 << U_force.torque, U_force.force;
                U_mat[i] = U_vec6;
                // D = Sᵀ·U = axis·U_torque (scalar stored as 1×1)
                S D = jt.axis.dot(U_force.torque);
                D_inv[i] = Eigen::Matrix<S,1,1>(S(1) / D);
                // u = τ - Sᵀ·p_A
                S pa_proj = jt.axis.dot(pA[bi].torque);
                u_vec[i] = Eigen::Matrix<S,1,1>(tau(model.v_idx[i]) - pa_proj);

                // I_A_bar = I_A - U · D⁻¹ · Uᵀ  (rank-1 update, 6×6)
                ArticulatedInertia<S> IA_bar;
                IA_bar.M = IA[bi].M - U_vec6 * D_inv[i](0,0) * U_vec6.transpose();

                // p_A_bar = p_A + I_A_bar · c + U · D⁻¹ · u
                Eigen::Matrix<S,6,1> c_vec;
                c_vec << c[i].angular, c[i].linear;
                Eigen::Matrix<S,6,1> pA_vec;
                pA_vec << pA[bi].torque, pA[bi].force;
                Eigen::Matrix<S,6,1> pA_bar_vec =
                    pA_vec + IA_bar.M * c_vec + U_vec6 * D_inv[i](0,0) * u_vec[i](0);
                Force6<S> pA_bar{pA_bar_vec.template head<3>(),
                                 pA_bar_vec.template tail<3>()};

                // Propagate to parent (skip world body — it has no dynamics).
                if (pi > 0) {
                    IA[pi] += X_up[i].transform_inertia(IA_bar);
                    pA[pi] += X_up[i].act_on_force(pA_bar);
                }

            } else {
                // FreeJoint: S = I_6, so U = I_A.M, D = I_A.M, I_A_bar = 0.
                Eigen::Matrix<S,6,1> tau_i = tau.segment(model.v_idx[i], 6);
                Eigen::Matrix<S,6,1> pA_vec;
                pA_vec << pA[bi].torque, pA[bi].force;

                U_mat[i] = IA[bi].M;             // 6×6
                D_inv[i] = IA[bi].M.inverse();   // 6×6
                u_vec[i] = tau_i - pA_vec;        // 6×1

                // I_A_bar = 0 (all DOFs free), so nothing propagates to parent.
                // (FreeJoint is typically at root; parent is world.)
            }
        }, model.joints[i]);
    }

    // ── Pass 3: compute accelerations root to tips ───────────────────────────
    Eigen::VectorX<S> ddq(model.nv);

    for (int i = 0; i < nj; ++i) {
        const int bi = i + 1;
        const int pi = model.parent[i];

        // Acceleration of child body before adding joint contribution.
        Motion6<S> a_before = X_up[i].inv_act_on_motion(a[pi]) + c[i];

        std::visit([&](const auto& jt) {
            using JT = std::decay_t<decltype(jt)>;

            if constexpr (std::is_same_v<JT, RevoluteJoint<S>>) {
                // ddq = D⁻¹ · (u - Uᵀ · a_before)
                // U is 6×1, a_before is Motion6 → cast to 6-vec.
                Eigen::Matrix<S,6,1> a_before_vec;
                a_before_vec << a_before.angular, a_before.linear;
                S ddq_i = D_inv[i](0,0) * (u_vec[i](0) - U_mat[i].col(0).dot(a_before_vec));
                ddq(model.v_idx[i]) = ddq_i;

                // Body acceleration: a_i = a_before + S * ddq_i
                Motion6<S> S_col{jt.axis, Eigen::Matrix<S,3,1>::Zero()};
                a[bi] = a_before + S_col * ddq_i;

            } else {
                // FreeJoint: ddq_6 = D⁻¹ · (u - I_A · a_before)
                Eigen::Matrix<S,6,1> a_before_vec;
                a_before_vec << a_before.angular, a_before.linear;
                Eigen::Matrix<S,6,1> ddq_6 =
                    D_inv[i] * (u_vec[i] - IA[bi].M * a_before_vec);
                ddq.segment(model.v_idx[i], 6) = ddq_6;

                // Body acceleration: a_i = a_before + S * ddq_6 = a_before + ddq_6 (S=I)
                a[bi].angular = a_before.angular + ddq_6.template head<3>();
                a[bi].linear  = a_before.linear  + ddq_6.template tail<3>();
            }
        }, model.joints[i]);
    }

    return ddq;
}

}  // namespace rapide
