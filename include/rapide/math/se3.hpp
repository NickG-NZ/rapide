#pragma once
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "rapide/math/spatial.hpp"

namespace rapide {

// SE3<S> — rigid transform: child frame relative to parent frame.
//
//   R : rotation (columns = child-frame axes in parent frame)
//   p : position of child-frame origin in parent frame
//
// Spatial actions are computed inline without materialising the 6×6 adjoint.
// The 6×6 adjoint is available via to_adjoint() for Jacobians / validation.
template<typename S>
struct SE3 {
    Eigen::Matrix<S,3,3> R;
    Eigen::Matrix<S,3,1> p;

    static SE3 Identity() {
        return {Eigen::Matrix<S,3,3>::Identity(), Eigen::Matrix<S,3,1>::Zero()};
    }

    // Compose: T_parent_child = T_parent_mid * T_mid_child
    SE3 operator*(const SE3& rhs) const {
        return {R * rhs.R, R * rhs.p + p};
    }

    // Inverse: T_child_parent
    SE3 inverse() const {
        Eigen::Matrix<S,3,3> Rt = R.transpose();
        return {Rt, -(Rt * p)};
    }

    // ── Spatial motion actions ───────────────────────────────────────────────
    // X (child→parent): ω_p = R·ω_c,  v_p = p×ω_p + R·v_c
    Motion6<S> act_on_motion(const Motion6<S>& m) const {
        Eigen::Matrix<S,3,1> omega = R * m.angular;
        return {omega, p.cross(omega) + R * m.linear};
    }

    // X⁻¹ (parent→child): ω_c = Rᵀ·ω_p,  v_c = Rᵀ·(v_p - p×ω_p)
    Motion6<S> inv_act_on_motion(const Motion6<S>& m) const {
        Eigen::Matrix<S,3,1> omega = R.transpose() * m.angular;
        return {omega, R.transpose() * (m.linear - p.cross(m.angular))};
    }

    // ── Spatial force actions ────────────────────────────────────────────────
    // X⁻ᵀ (child→parent, dual of X⁻¹ for motion):
    //   f_p = R·f_c,  τ_p = R·τ_c - p×f_p
    Force6<S> act_on_force(const Force6<S>& f) const {
        Eigen::Matrix<S,3,1> flin = R * f.force;
        return {R * f.torque - p.cross(flin), flin};
    }

    // Xᵀ (parent→child, dual of X for motion):
    //   f_c = Rᵀ·f_p,  τ_c = Rᵀ·(τ_p + p×f_p)
    Force6<S> inv_act_on_force(const Force6<S>& f) const {
        Eigen::Matrix<S,3,1> flin = R.transpose() * f.force;
        return {R.transpose() * (f.torque + p.cross(f.force)), flin};
    }

    // ── Inertia transform ────────────────────────────────────────────────────
    // Transform ArticulatedInertia from child frame to parent frame:
    //   I_parent = X⁻ᵀ · I_child · X⁻¹
    // Uses the materialised 6×6 adjoint since ArticulatedInertia is already 6×6.
    ArticulatedInertia<S> transform_inertia(const ArticulatedInertia<S>& I_child) const {
        Eigen::Matrix<S,6,6> Xinv = to_inv_adjoint();
        return {Xinv.transpose() * I_child.M * Xinv};
    }

    // ── Materialise 6×6 (escape hatch) ──────────────────────────────────────
    // Adjoint of X (motion child→parent):
    //   [ R,    0  ]
    //   [ p×R,  R  ]
    Eigen::Matrix<S,6,6> to_adjoint() const {
        Eigen::Matrix<S,6,6> X;
        X.setZero();
        X.template topLeftCorner<3,3>()     = R;
        X.template bottomLeftCorner<3,3>()  = skew(p) * R;
        X.template bottomRightCorner<3,3>() = R;
        return X;
    }

    // Adjoint of X⁻¹ (motion parent→child):
    //   [ Rᵀ,       0  ]
    //   [ -Rᵀ·p×,  Rᵀ  ]
    Eigen::Matrix<S,6,6> to_inv_adjoint() const {
        Eigen::Matrix<S,6,6> X;
        X.setZero();
        Eigen::Matrix<S,3,3> Rt = R.transpose();
        X.template topLeftCorner<3,3>()     =  Rt;
        X.template bottomLeftCorner<3,3>()  = -Rt * skew(p);
        X.template bottomRightCorner<3,3>() =  Rt;
        return X;
    }

    static Eigen::Matrix<S,3,3> skew(const Eigen::Matrix<S,3,1>& v) {
        Eigen::Matrix<S,3,3> m;
        m <<  S(0), -v(2),  v(1),
              v(2),  S(0), -v(0),
             -v(1),  v(0),  S(0);
        return m;
    }
};

// Construct SE3 from a unit quaternion + position.
template<typename S>
SE3<S> se3_from_quat_pos(const Eigen::Quaternion<S>& q,
                          const Eigen::Matrix<S,3,1>& p) {
    return {q.toRotationMatrix(), p};
}

}  // namespace rapide
