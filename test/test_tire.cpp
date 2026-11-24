#include "rapide/contact/tire.hpp"
#include "rapide/math/se3.hpp"
#include <catch2/catch_test_macros.hpp>

using namespace rapide;
using S = double;
using V3 = Eigen::Matrix<S, 3, 1>;
using M3 = Eigen::Matrix<S, 3, 3>;

static constexpr S eps = 1e-9;

TEST_CASE("TireContact: no contact when wheel is above ground") {
  TireContact<S> tire{0.3, 1e5, 1e3, 0.8};
  SE3<S> X{M3::Identity(), V3{0, 0, 0.5}};
  Twist<S> v = Twist<S>::Zero();
  Wrench<S> f = tire.compute(X, v);
  REQUIRE(f.force.norm() < eps);
  REQUIRE(f.torque.norm() < eps);
}

TEST_CASE("TireContact: normal force increases with penetration") {
  TireContact<S> tire{0.3, 1e5, 1e3, 0.8};
  Twist<S> v = Twist<S>::Zero();

  SE3<S> X1{M3::Identity(), V3{0, 0, 0.2}}; // penetration = 0.1
  SE3<S> X2{M3::Identity(), V3{0, 0, 0.1}}; // penetration = 0.2

  Wrench<S> f1 = tire.compute(X1, v);
  Wrench<S> f2 = tire.compute(X2, v);

  REQUIRE(f2.force(2) > f1.force(2));
  REQUIRE(f1.force(2) > 0.0);
}

TEST_CASE("TireContact: friction opposes forward slip velocity") {
  TireContact<S> tire{0.3, 1e5, 1e3, 0.8};
  SE3<S> X{M3::Identity(),
           V3{0, 0, 0.25}}; // penetrating (radius=0.3, centre at 0.25)
  Twist<S> v;
  v.angular = V3::Zero();
  v.linear = V3{1.0, 0.0, 0.0}; // moving forward in body frame

  Wrench<S> f = tire.compute(X, v);
  REQUIRE(f.force(0) < 0.0);           // friction opposes +X motion
  REQUIRE(std::abs(f.force(1)) < eps); // no Y slip → no Y friction
}

TEST_CASE(
    "TireContact: normal force clamped to zero when spring-damper pulls") {
  TireContact<S> tire{0.3, 100.0, 1e6, 0.8}; // high damping, low stiffness
  SE3<S> X{M3::Identity(), V3{0, 0, 0.29}};  // tiny penetration
  Twist<S> v;
  v.angular = V3::Zero();
  v.linear = V3{0.0, 0.0, 10.0}; // moving up fast

  Wrench<S> f = tire.compute(X, v);
  REQUIRE(f.force(2) >= 0.0); // clamped — no suction
}
