#pragma once
#include <Eigen/Core>

namespace rapide {

// Spatial motion vector [angular (ω); linear (v)]  — Featherstone convention
template<typename S>
struct Motion6 {
    Eigen::Matrix<S,3,1> angular;
    Eigen::Matrix<S,3,1> linear;

    static Motion6 Zero() {
        return {Eigen::Matrix<S,3,1>::Zero(), Eigen::Matrix<S,3,1>::Zero()};
    }

    Motion6 operator+(const Motion6& o) const { return {angular + o.angular, linear + o.linear}; }
    Motion6 operator-(const Motion6& o) const { return {angular - o.angular, linear - o.linear}; }
    Motion6& operator+=(const Motion6& o) { angular += o.angular; linear += o.linear; return *this; }
    Motion6 operator*(S s) const { return {angular * s, linear * s}; }

    // Spatial cross product for motion × motion:  [this × m]
    // (v₁ ×) v₂ = [ω₁×ω₂; ω₁×v₂ + v₁×ω₂]
    Motion6 cross(const Motion6& m) const {
        return {
            angular.cross(m.angular),
            angular.cross(m.linear) + linear.cross(m.angular)
        };
    }
};

// Spatial force vector [torque (τ); force (f)]
template<typename S>
struct Force6 {
    Eigen::Matrix<S,3,1> torque;
    Eigen::Matrix<S,3,1> force;

    static Force6 Zero() {
        return {Eigen::Matrix<S,3,1>::Zero(), Eigen::Matrix<S,3,1>::Zero()};
    }

    Force6 operator+(const Force6& o) const { return {torque + o.torque, force + o.force}; }
    Force6 operator-(const Force6& o) const { return {torque - o.torque, force - o.force}; }
    Force6& operator+=(const Force6& o) { torque += o.torque; force += o.force; return *this; }
    Force6 operator*(S s) const { return {torque * s, force * s}; }

    S dot(const Motion6<S>& m) const { return torque.dot(m.angular) + force.dot(m.linear); }
};

// Dual cross product: motion acting on force   (v ×*) f
// Used for ABA Coriolis bias:  p += v ×* I*v
template<typename S>
Force6<S> cross_motion_force(const Motion6<S>& v, const Force6<S>& f) {
    return {
        v.angular.cross(f.torque) + v.linear.cross(f.force),
        v.angular.cross(f.force)
    };
}

// ─── Rigid-body spatial inertia ──────────────────────────────────────────────
// Stored compactly as (mass, com, I_c). The 6×6 matrix is:
//   [ I_c - m·c×², -m·c× ]
//   [   m·c×,       m·I  ]
// where c× = skew(c). Never stored in that form; computed on demand.
template<typename S>
struct SpatialInertia {
    S                      mass;
    Eigen::Matrix<S,3,1>   com;   // centre of mass in body frame
    Eigen::Matrix<S,3,3>   I_c;   // rotational inertia about CoM

    // I * m — apply inertia to a motion vector, return a force vector.
    // f_lin = m·(c×ω + v)
    // f_ang = I_c·ω - c × f_lin
    Force6<S> apply(const Motion6<S>& m) const {
        Eigen::Matrix<S,3,1> f_lin = mass * (com.cross(m.angular) + m.linear);
        Eigen::Matrix<S,3,1> f_ang = I_c * m.angular - com.cross(f_lin);
        return {f_ang, f_lin};
    }

    // Materialise the 6×6 for validation or explicit Jacobians.
    Eigen::Matrix<S,6,6> to_matrix() const {
        const Eigen::Matrix<S,3,3> cx = skew(com);
        Eigen::Matrix<S,6,6> M;
        // Top-left:  I_c - m·cx²  = I_c + m·cx·cx^T
        M.template topLeftCorner<3,3>()     = I_c - mass * cx * cx;
        // Top-right: -m·cx
        M.template topRightCorner<3,3>()    = -mass * cx;
        // Bottom-left:  m·cx
        M.template bottomLeftCorner<3,3>()  =  mass * cx;
        // Bottom-right: m·I
        M.template bottomRightCorner<3,3>() =  mass * Eigen::Matrix<S,3,3>::Identity();
        return M;
    }

    static Eigen::Matrix<S,3,3> skew(const Eigen::Matrix<S,3,1>& v) {
        Eigen::Matrix<S,3,3> m;
        m <<  S(0), -v(2),  v(1),
              v(2),  S(0), -v(0),
             -v(1),  v(0),  S(0);
        return m;
    }
};

// ─── Articulated-body inertia ─────────────────────────────────────────────────
// Used internally by ABA. Cannot be stored compactly in general; full 6×6.
template<typename S>
struct ArticulatedInertia {
    Eigen::Matrix<S,6,6> M;

    static ArticulatedInertia Zero() {
        ArticulatedInertia a;
        a.M.setZero();
        return a;
    }

    ArticulatedInertia operator+(const ArticulatedInertia& o) const { return {M + o.M}; }
    ArticulatedInertia& operator+=(const ArticulatedInertia& o) { M += o.M; return *this; }

    Force6<S> apply(const Motion6<S>& m) const {
        Eigen::Matrix<S,6,1> mv;
        mv << m.angular, m.linear;
        Eigen::Matrix<S,6,1> fv = M * mv;
        return {fv.template head<3>(), fv.template tail<3>()};
    }

    static ArticulatedInertia from_rigid(const SpatialInertia<S>& I) {
        return {I.to_matrix()};
    }
};

}  // namespace rapide
