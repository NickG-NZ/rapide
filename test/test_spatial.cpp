#include "rapide/math/spatial.hpp"
#include <catch2/catch_test_macros.hpp>

using namespace rapide;
using S = double;
using V3 = Eigen::Matrix<S, 3, 1>;
using M3 = Eigen::Matrix<S, 3, 3>;

static constexpr S eps = 1e-10;

TEST_CASE("SpatialInertia::apply matches 6x6 matrix multiplication") {
  // Construct a known inertia.
  S mass = 3.0;
  V3 com{0.1, -0.2, 0.3};
  M3 I_c;
  I_c << 2.0, 0.1, 0.0, 0.1, 1.5, 0.0, 0.0, 0.0, 2.5;

  SpatialInertia<S> I{mass, com, I_c};
  Twist<S> m{{1.0, 0.5, -0.3}, {0.2, -0.1, 0.4}};

  Wrench<S> f_method = I.apply(m);

  // Compare against full 6×6 multiply.
  Eigen::Matrix<S, 6, 6> M6 = I.to_matrix();
  Eigen::Matrix<S, 6, 1> mv;
  mv << m.angular, m.linear;
  Eigen::Matrix<S, 6, 1> fv = M6 * mv;

  Eigen::Matrix<S, 6, 1> f_method_vec;
  f_method_vec << f_method.torque, f_method.force;
  REQUIRE((f_method_vec - fv).norm() < eps);
}

TEST_CASE("SpatialInertia::to_matrix is symmetric") {
  SpatialInertia<S> I{2.0, {0.1, 0.0, 0.2}, M3::Identity() * 0.5};
  Eigen::Matrix<S, 6, 6> M = I.to_matrix();
  REQUIRE((M - M.transpose()).norm() < eps);
}

TEST_CASE("SpatialInertia kinetic energy is non-negative") {
  SpatialInertia<S> I{5.0, {0.0, 0.0, 0.1}, M3::Identity()};
  Twist<S> m{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}};
  // KE = 0.5 * m^T * I * m
  Wrench<S> Im = I.apply(m);
  S ke = 0.5 * Im.dot(m);
  REQUIRE(ke > 0.0);
}

TEST_CASE("cross_motion_force matches 6x6 computation") {
  Twist<S> v{{1.0, 0.0, 0.5}, {0.0, 1.0, 0.0}};
  Wrench<S> f{{0.5, -0.3, 0.1}, {1.0, 0.0, -0.2}};

  Wrench<S> result = cross_motion_force(v, f);

  // Build the 6×6 spatial cross operator for motion and apply to force.
  // (v ×*) = -（v ×）^T in Featherstone's notation.
  // Manual check: {ω×τ + v×f; ω×f}
  V3 expected_torque = v.angular.cross(f.torque) + v.linear.cross(f.force);
  V3 expected_force = v.angular.cross(f.force);

  REQUIRE((result.torque - expected_torque).norm() < eps);
  REQUIRE((result.force - expected_force).norm() < eps);
}

TEST_CASE("ArticulatedInertia::from_rigid recovers 6x6") {
  SpatialInertia<S> rigid{4.0, {0.1, 0.2, 0.0}, M3::Identity() * 2.0};
  ArticulatedInertia<S> art = ArticulatedInertia<S>::from_rigid(rigid);
  REQUIRE((art.M - rigid.to_matrix()).norm() < eps);
}

TEST_CASE("Offset point mass has physical linear and angular momentum") {
  SpatialInertia<double> inertia{1., {1, 0, 0}, M3::Zero()};
  Twist<double> velocity{{0, 0, 1}, {0, 1, 0}};
  const auto momentum = inertia.apply(velocity);
  REQUIRE((momentum.force - V3{0, 2, 0}).norm() < eps);
  REQUIRE((momentum.torque - V3{0, 0, 2}).norm() < eps);
  REQUIRE(std::abs(.5 * momentum.dot(velocity) - 2.) < eps);
}
