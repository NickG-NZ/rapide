#pragma once
#include "rapide/aba/tree_data.hpp"
#include "rapide/math/se3.hpp"
#include "rapide/math/spatial.hpp"
#include <Eigen/Core>
#include <Eigen/LU>
#include <vector>

namespace rapide::aba {

/** Return generalized accelerations using Featherstone's articulated-body
 *  algorithm without forming a global mass matrix. A forward pass computes
 *  transforms, velocities, and bias forces; a backward pass eliminates joint
 *  accelerations and accumulates articulated inertias; a final forward pass
 *  recovers accelerations from the root to the leaves.
 *
 *  tau contains generalized joint efforts. f_ext contains world-expressed
 *  spatial wrenches about the world origin, indexed by body (including world).
 *  Gravity enters through the base acceleration and must not also be in f_ext.
 *  Floating joints are supported at roots; fixed joints propagate inertia
 *  without adding generalized coordinates. */
template <typename S>
Eigen::VectorX<S> forward_dynamics(
    const TreeData<S> &model, const JointState<S> &state,
    const Eigen::VectorX<S> &tau, const std::vector<Wrench<S>> &f_ext,
    const Eigen::Matrix<S, 3, 1> &gravity = Eigen::Matrix<S, 3, 1>(S(0), S(0),
                                                                   S(-9.81))) {
  const int nb = model.num_bodies();
  const int nj = model.num_joints();

  // ── Workspace ─────────────────────────────────────────────────────────────
  std::vector<SE3<S>> X_up(nj);
  std::vector<SE3<S>> X_world(nb, SE3<S>::Identity());
  std::vector<Twist<S>> v(nb, Twist<S>::Zero());
  std::vector<Twist<S>> c(nj, Twist<S>::Zero());
  std::vector<ArticulatedInertia<S>> IA(nb, ArticulatedInertia<S>::Zero());
  std::vector<Wrench<S>> pA(nb, Wrench<S>::Zero());
  std::vector<Twist<S>> a(nb, Twist<S>::Zero());

  // Per-joint data carried from pass 2 → pass 3.
  std::vector<Eigen::MatrixX<S>> D_inv(nj); // nv_i × nv_i
  std::vector<Eigen::VectorX<S>> u_vec(nj); // nv_i
  std::vector<Eigen::MatrixX<S>> U_mat(nj); // 6 × nv_i

  // Gravity trick (RBDA §9.3): setting a[world] = -g makes gravity appear
  // automatically as an inertial bias in every body.
  a[0] = {Eigen::Matrix<S, 3, 1>::Zero(), -gravity};

  // ── Pass 1: kinematics ────────────────────────────────────────────────────
  for (int i = 0; i < nj; ++i) {
    const int bi = i + 1;
    const int pi = model.parent[i];

    SE3<S> X_J;
    Twist<S> vJ;
    std::visit(
        [&](const auto &jt) {
          using JT = std::decay_t<decltype(jt)>;
          if constexpr (std::is_same_v<JT, FreeJoint<S>>) {
            Eigen::Matrix<S, 7, 1> qi = model.q_i(state.q, i);
            Eigen::Matrix<S, 6, 1> vi = model.v_i(state.v, i);
            X_J = jt.X_J(qi);
            vJ = jt.vJ(vi);
          } else if constexpr (std::is_same_v<JT, FixedJoint<S>>) {
            X_J = jt.X_J();
            vJ = Twist<S>::Zero();
          } else {
            X_J = jt.X_J(state.q(model.q_idx[i]));
            vJ = jt.vJ(state.v(model.v_idx[i]));
          }
        },
        model.joints[i]);

    X_up[i] = model.X_tree[i] * X_J;
    X_world[bi] = X_world[pi] * X_up[i];
    v[bi] = X_up[i].inv_act_on_motion(v[pi]) + vJ;
    c[i] = v[bi].cross(vJ);

    IA[bi] = ArticulatedInertia<S>::from_rigid(model.bodies[bi].inertia);
    Wrench<S> fext_body = X_world[bi].inv_act_on_force(f_ext[bi]);
    pA[bi] = cross_motion_force(v[bi], model.bodies[bi].inertia.apply(v[bi])) -
             fext_body;
  }

  // ── Pass 2: tip → root, accumulate articulated inertia & bias ────────────
  for (int i = nj - 1; i >= 0; --i) {
    const int bi = i + 1;
    const int pi = model.parent[i];

    std::visit(
        [&](const auto &jt) {
          using JT = std::decay_t<decltype(jt)>;

          if constexpr (std::is_same_v<JT, FreeJoint<S>>) {
            // 6-DOF: S = I₆, so U = IA.M, D = IA.M, I_A_bar = 0.
            Eigen::Matrix<S, 6, 1> pA_vec;
            pA_vec << pA[bi].torque, pA[bi].force;
            U_mat[i] = IA[bi].M;
            D_inv[i] = IA[bi].M.inverse();
            u_vec[i] = tau.segment(model.v_idx[i], 6) - pA_vec;
            // I_A_bar = 0 → nothing propagates to parent (parent = world).

          } else if constexpr (std::is_same_v<JT, FixedJoint<S>>) {
            if (pi > 0) {
              IA[pi] += X_up[i].transform_inertia(IA[bi]);
              pA[pi] += X_up[i].act_on_force(pA[bi] + IA[bi].apply(c[i]));
            }
          } else {
            // Generic 1-DOF: works for RevoluteJoint and PrismaticJoint alike.
            // Both provide S_times_v(1) and S_transpose_force().
            Twist<S> S_col = jt.S_times_v(S(1));
            Wrench<S> U_force = IA[bi].apply(S_col);
            Eigen::Matrix<S, 6, 1> U_vec6;
            U_vec6 << U_force.torque, U_force.force;
            U_mat[i] = U_vec6;

            S D = jt.S_transpose_force(U_force);
            D_inv[i] = Eigen::Matrix<S, 1, 1>(S(1) / D);

            S pa_proj = jt.S_transpose_force(pA[bi]);
            u_vec[i] = Eigen::Matrix<S, 1, 1>(tau(model.v_idx[i]) - pa_proj);

            // I_A_bar = I_A − U·D⁻¹·Uᵀ  (rank-1 downdate)
            ArticulatedInertia<S> IA_bar;
            IA_bar.M = IA[bi].M - U_vec6 * D_inv[i](0, 0) * U_vec6.transpose();

            // p_A_bar = p_A + I_A_bar·c + U·D⁻¹·u
            Eigen::Matrix<S, 6, 1> c_vec;
            c_vec << c[i].angular, c[i].linear;
            Eigen::Matrix<S, 6, 1> pA_vec;
            pA_vec << pA[bi].torque, pA[bi].force;
            Eigen::Matrix<S, 6, 1> pA_bar_vec =
                pA_vec + IA_bar.M * c_vec +
                U_vec6 * D_inv[i](0, 0) * u_vec[i](0);
            Wrench<S> pA_bar{pA_bar_vec.template head<3>(),
                             pA_bar_vec.template tail<3>()};

            if (pi > 0) {
              IA[pi] += X_up[i].transform_inertia(IA_bar);
              pA[pi] += X_up[i].act_on_force(pA_bar);
            }
          }
        },
        model.joints[i]);
  }

  // ── Pass 3: root → tip, resolve accelerations ─────────────────────────────
  Eigen::VectorX<S> ddq(model.nv);

  for (int i = 0; i < nj; ++i) {
    const int bi = i + 1;
    const int pi = model.parent[i];

    Twist<S> a_before = X_up[i].inv_act_on_motion(a[pi]) + c[i];

    std::visit(
        [&](const auto &jt) {
          using JT = std::decay_t<decltype(jt)>;

          if constexpr (std::is_same_v<JT, FreeJoint<S>>) {
            Eigen::Matrix<S, 6, 1> a_bv;
            a_bv << a_before.angular, a_before.linear;
            Eigen::Matrix<S, 6, 1> ddq_6 =
                D_inv[i] * (u_vec[i] - IA[bi].M * a_bv);
            ddq.segment(model.v_idx[i], 6) = ddq_6;
            a[bi].angular = a_before.angular + ddq_6.template head<3>();
            a[bi].linear = a_before.linear + ddq_6.template tail<3>();

          } else if constexpr (std::is_same_v<JT, FixedJoint<S>>) {
            a[bi] = a_before;
          } else {
            // Generic 1-DOF
            Eigen::Matrix<S, 6, 1> a_bv;
            a_bv << a_before.angular, a_before.linear;
            S ddq_i =
                D_inv[i](0, 0) * (u_vec[i](0) - U_mat[i].col(0).dot(a_bv));
            ddq(model.v_idx[i]) = ddq_i;
            a[bi] = a_before + jt.S_times_v(ddq_i);
          }
        },
        model.joints[i]);
  }

  return ddq;
}

} // namespace rapide::aba
