#include <catch2/catch_test_macros.hpp>
#include "rapide/math/se3.hpp"

using namespace rapide;
using S = double;
using V3 = Eigen::Matrix<S,3,1>;
using M3 = Eigen::Matrix<S,3,3>;

static constexpr S eps = 1e-10;

TEST_CASE("SE3 identity is a no-op on motion vectors") {
    SE3<S> I = SE3<S>::Identity();
    Motion6<S> m{{1, 2, 3}, {4, 5, 6}};
    auto m2 = I.act_on_motion(m);
    REQUIRE((m2.angular - m.angular).norm() < eps);
    REQUIRE((m2.linear  - m.linear ).norm() < eps);
}

TEST_CASE("SE3 compose with inverse yields identity") {
    Eigen::AngleAxis<S> aa(0.4, V3(1,1,0).normalized());
    SE3<S> T{aa.toRotationMatrix(), V3(1, 2, 3)};
    SE3<S> prod = T * T.inverse();
    REQUIRE((prod.R - M3::Identity()).norm() < eps);
    REQUIRE(prod.p.norm() < eps);
}

TEST_CASE("SE3 act_on_motion and inv_act_on_motion are inverses") {
    Eigen::AngleAxis<S> aa(1.1, V3(0,1,0).normalized());
    SE3<S> T{aa.toRotationMatrix(), V3(1, -1, 2)};

    Motion6<S> m{{0.5, -0.3, 0.2}, {1.0, 0.0, -0.5}};
    Motion6<S> m_back = T.inv_act_on_motion(T.act_on_motion(m));

    REQUIRE((m_back.angular - m.angular).norm() < eps);
    REQUIRE((m_back.linear  - m.linear ).norm() < eps);
}

TEST_CASE("SE3 to_adjoint matches act_on_motion element-wise") {
    Eigen::AngleAxis<S> aa(0.9, V3(1,2,3).normalized());
    SE3<S> T{aa.toRotationMatrix(), V3(-1, 2, 0.5)};

    Motion6<S> m{{0.1, -0.2, 0.3}, {-1, 0.5, 0}};
    Motion6<S> m_method = T.act_on_motion(m);

    Eigen::Matrix<S,6,1> mv; mv << m.angular, m.linear;
    Eigen::Matrix<S,6,1> mv2 = T.to_adjoint() * mv;
    Eigen::Matrix<S,6,1> m_method_vec; m_method_vec << m_method.angular, m_method.linear;
    REQUIRE((m_method_vec - mv2).norm() < eps);
}

TEST_CASE("SE3 to_adjoint and to_inv_adjoint are matrix inverses") {
    Eigen::AngleAxis<S> aa(0.7, V3(1,0,0).normalized());
    SE3<S> T{aa.toRotationMatrix(), V3(0.5, 1.5, -1.0)};
    auto X    = T.to_adjoint();
    auto Xinv = T.to_inv_adjoint();
    REQUIRE((X * Xinv - Eigen::Matrix<S,6,6>::Identity()).norm() < eps);
}

TEST_CASE("SE3 power invariance: force-motion inner product preserved across frames") {
    Eigen::AngleAxis<S> aa(0.7, V3(1,0,0).normalized());
    SE3<S> T{aa.toRotationMatrix(), V3(0.5, 1.5, -1.0)};

    Motion6<S> m{{1, 0, 0}, {0, 1, 0}};
    Force6<S>  f{{0, 1, 0}, {0, 0, 1}};

    // Power in child frame
    S power_child = f.dot(m);

    // Transform to parent frame: X^{-T} for force, X for motion (X = to_adjoint)
    auto X    = T.to_adjoint();
    auto Xinv = T.to_inv_adjoint();
    Eigen::Matrix<S,6,1> fv; fv << f.torque, f.force;
    Eigen::Matrix<S,6,1> mv; mv << m.angular, m.linear;

    S power_parent = (Xinv.transpose() * fv).dot(X * mv);
    REQUIRE(std::abs(power_child - power_parent) < eps);
}
