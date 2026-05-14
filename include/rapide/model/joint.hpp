#pragma once
#include <variant>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "rapide/math/se3.hpp"
#include "rapide/math/spatial.hpp"

namespace rapide {

// ─── RevoluteJoint ───────────────────────────────────────────────────────────
// Single-axis rotation. nq = nv = 1.
// axis: unit vector in joint predecessor frame.
template<typename S>
struct RevoluteJoint {
    Eigen::Matrix<S,3,1> axis;  // rotation axis in predecessor frame

    static constexpr int NQ = 1;
    static constexpr int NV = 1;

    // Transform across joint: rotation by angle q about axis.
    SE3<S> X_J(S q) const {
        Eigen::AngleAxis<S> aa(q, axis);
        return {aa.toRotationMatrix(), Eigen::Matrix<S,3,1>::Zero()};
    }

    // Motion subspace column (the 6×1 joint motion axis in joint frame).
    // For a revolute joint this is constant: S = [axis; 0].
    Motion6<S> S_times_v(S dq) const {
        return {axis * dq, Eigen::Matrix<S,3,1>::Zero()};
    }

    // Joint velocity contribution (vJ = S * dq).
    Motion6<S> vJ(S dq) const { return S_times_v(dq); }

    // S^T * f: project force onto the joint axis (scalar output for 1-DOF).
    S S_transpose_force(const Force6<S>& f) const {
        return axis.dot(f.torque);
    }

    // Integrate configuration: q_new = q + dq * dt (trivial for revolute).
    S integrate_q(S q, S dq, S dt) const { return q + dq * dt; }
};

// ─── FreeJoint ────────────────────────────────────────────────────────────────
// 6-DOF floating base. nq = 7 (quaternion xyzw + position xyz), nv = 6.
// Velocity is body-frame spatial velocity [ω; v].
template<typename S>
struct FreeJoint {
    static constexpr int NQ = 7;
    static constexpr int NV = 6;

    // q_i: [qx, qy, qz, qw, px, py, pz]
    SE3<S> X_J(const Eigen::Ref<const Eigen::Matrix<S,7,1>>& q) const {
        Eigen::Quaternion<S> quat(q(3), q(0), q(1), q(2));  // (w, x, y, z)
        quat.normalize();
        return {quat.toRotationMatrix(), q.template tail<3>()};
    }

    // vJ = [ω_body; v_body] — the full 6D body-frame velocity.
    Motion6<S> vJ(const Eigen::Ref<const Eigen::Matrix<S,6,1>>& dq) const {
        return {dq.template head<3>(), dq.template tail<3>()};
    }

    // S^T * f: all 6 components (returns the full generalised force vector).
    Eigen::Matrix<S,6,1> S_transpose_force(const Force6<S>& f) const {
        Eigen::Matrix<S,6,1> tau;
        tau << f.torque, f.force;
        return tau;
    }

    // Integrate configuration using first-order quaternion kinematics:
    //   q_new = q ⊗ exp(ω * dt / 2),  then renormalise
    //   p_new = p + R * v_body * dt
    Eigen::Matrix<S,7,1> integrate_q(
        const Eigen::Ref<const Eigen::Matrix<S,7,1>>& q,
        const Eigen::Ref<const Eigen::Matrix<S,6,1>>& dq,
        S dt) const
    {
        Eigen::Quaternion<S> quat(q(3), q(0), q(1), q(2));
        Eigen::Matrix<S,3,1> omega = dq.template head<3>();
        Eigen::Matrix<S,3,1> v_body = dq.template tail<3>();

        // First-order: delta_q = [0, ω*dt/2] in quaternion notation
        S half_dt = dt / S(2);
        Eigen::Quaternion<S> delta(S(1), omega(0)*half_dt, omega(1)*half_dt, omega(2)*half_dt);
        Eigen::Quaternion<S> quat_new = (quat * delta).normalized();

        Eigen::Matrix<S,3,1> p = q.template tail<3>();
        Eigen::Matrix<S,3,1> p_new = p + quat.toRotationMatrix() * v_body * dt;

        Eigen::Matrix<S,7,1> q_new;
        q_new << quat_new.x(), quat_new.y(), quat_new.z(), quat_new.w(), p_new;
        return q_new;
    }
};

// ─── Joint variant ────────────────────────────────────────────────────────────
template<typename S>
using Joint = std::variant<FreeJoint<S>, RevoluteJoint<S>>;

// Accessors that work via std::visit.
template<typename S>
int joint_nq(const Joint<S>& j) {
    return std::visit([](const auto& jt) { return jt.NQ; }, j);
}

template<typename S>
int joint_nv(const Joint<S>& j) {
    return std::visit([](const auto& jt) { return jt.NV; }, j);
}

}  // namespace rapide
