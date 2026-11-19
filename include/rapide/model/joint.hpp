#pragma once
#include <variant>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "rapide/math/se3.hpp"
#include "rapide/math/spatial.hpp"

namespace rapide {

// ─── RevoluteJoint ───────────────────────────────────────────────────────────
// Single-axis rotation. nq = nv = 1.
template<typename S>
struct RevoluteJoint {
    Eigen::Matrix<S,3,1> axis;  // unit rotation axis in predecessor frame

    static constexpr int NQ = 1;
    static constexpr int NV = 1;

    SE3<S> X_J(S q) const {
        Eigen::AngleAxis<S> aa(q, axis);
        return {aa.toRotationMatrix(), Eigen::Matrix<S,3,1>::Zero()};
    }

    // Motion subspace: S_col = [axis; 0]
    Twist<S> S_times_v(S dq) const {
        return {axis * dq, Eigen::Matrix<S,3,1>::Zero()};
    }
    Twist<S> vJ(S dq) const { return S_times_v(dq); }

    // S^T * f = axis · f.torque
    S S_transpose_force(const Wrench<S>& f) const { return axis.dot(f.torque); }

    S integrate_q(S q, S dq, S dt) const { return q + dq * dt; }
};

// ─── PrismaticJoint ──────────────────────────────────────────────────────────
// Single-axis translation. nq = nv = 1.
template<typename S>
struct PrismaticJoint {
    Eigen::Matrix<S,3,1> axis;  // unit translation axis in predecessor frame

    static constexpr int NQ = 1;
    static constexpr int NV = 1;

    SE3<S> X_J(S q) const {
        return {Eigen::Matrix<S,3,3>::Identity(), axis * q};
    }

    // Motion subspace: S_col = [0; axis]
    Twist<S> S_times_v(S dq) const {
        return {Eigen::Matrix<S,3,1>::Zero(), axis * dq};
    }
    Twist<S> vJ(S dq) const { return S_times_v(dq); }

    // S^T * f = axis · f.force   (translational — NOT torque)
    S S_transpose_force(const Wrench<S>& f) const { return axis.dot(f.force); }

    S integrate_q(S q, S dq, S dt) const { return q + dq * dt; }
};

// ─── FreeJoint ────────────────────────────────────────────────────────────────
// 6-DOF floating base. nq = 7 (quaternion xyzw + position xyz), nv = 6.
// Velocity is body-frame spatial velocity [ω; v].
template<typename S>
struct FreeJoint {
    static constexpr int NQ = 7;
    static constexpr int NV = 6;

    SE3<S> X_J(const Eigen::Ref<const Eigen::Matrix<S,7,1>>& q) const {
        Eigen::Quaternion<S> quat(q(3), q(0), q(1), q(2));  // (w, x, y, z)
        quat.normalize();
        return {quat.toRotationMatrix(), q.template tail<3>()};
    }

    Twist<S> vJ(const Eigen::Ref<const Eigen::Matrix<S,6,1>>& dq) const {
        return {dq.template head<3>(), dq.template tail<3>()};
    }

    Eigen::Matrix<S,6,1> S_transpose_force(const Wrench<S>& f) const {
        Eigen::Matrix<S,6,1> tau;
        tau << f.torque, f.force;
        return tau;
    }

    // Integrate quaternion via first-order kinematics, then renormalise.
    Eigen::Matrix<S,7,1> integrate_q(
        const Eigen::Ref<const Eigen::Matrix<S,7,1>>& q,
        const Eigen::Ref<const Eigen::Matrix<S,6,1>>& dq,
        S dt) const
    {
        Eigen::Quaternion<S> quat(q(3), q(0), q(1), q(2));
        Eigen::Matrix<S,3,1> omega  = dq.template head<3>();
        Eigen::Matrix<S,3,1> v_body = dq.template tail<3>();

        S half_dt = dt / S(2);
        Eigen::Quaternion<S> delta(S(1),
            omega(0)*half_dt, omega(1)*half_dt, omega(2)*half_dt);
        Eigen::Quaternion<S> quat_new = (quat * delta).normalized();

        Eigen::Matrix<S,3,1> p_new =
            q.template tail<3>() + quat.toRotationMatrix() * v_body * dt;

        Eigen::Matrix<S,7,1> q_new;
        q_new << quat_new.x(), quat_new.y(), quat_new.z(), quat_new.w(), p_new;
        return q_new;
    }
};

// ─── Joint variant ────────────────────────────────────────────────────────────
template<typename S>
using Joint = std::variant<FreeJoint<S>, RevoluteJoint<S>, PrismaticJoint<S>>;

template<typename S>
int joint_nq(const Joint<S>& j) {
    return std::visit([](const auto& jt) { return jt.NQ; }, j);
}

template<typename S>
int joint_nv(const Joint<S>& j) {
    return std::visit([](const auto& jt) { return jt.NV; }, j);
}

}  // namespace rapide
