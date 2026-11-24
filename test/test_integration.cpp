#include "rapide/integration/euler.hpp"
#include "rapide/integration/joint_configuration.hpp"
#include "rapide/simulator.hpp"
#include <catch2/catch_test_macros.hpp>
#include <numbers>

TEST_CASE("Shared orientation update applies exact body-axis delta rotations") {
  using namespace rapide;
  const Eigen::Quaterniond initial(Eigen::AngleAxisd(.7, Vec3::UnitX()));
  const Vec3 omega{0, 0, std::numbers::pi};
  const auto expected =
      initial *
      Eigen::Quaterniond(Eigen::AngleAxisd(std::numbers::pi, Vec3::UnitZ()));
  const auto next = integration::integrate_orientation(initial, omega, 1.);
  REQUIRE(next.angularDistance(expected) < 1e-12);
  REQUIRE(integration::integrate_orientation(initial, Vec3::Zero(), 1.)
              .angularDistance(initial) < 1e-12);
  const Vec3 tiny{1e-14, -2e-14, 3e-14};
  REQUIRE((rotation_increment(tiny).vec() - .5 * tiny).norm() < 1e-25);
  auto repeated = initial;
  for (int i = 0; i < 10000; ++i) {
    repeated = integration::integrate_orientation(repeated, omega, .0001);
  }
  REQUIRE(std::abs(repeated.norm() - 1) < 1e-12);
  REQUIRE(repeated.angularDistance(expected) < 1e-11);
  const Eigen::Vector3f float_omega{0, 0, 1};
  const auto f = integration::integrate_orientation(
      Eigen::Quaternionf::Identity(), float_omega, .5f);
  REQUIRE(f.angularDistance(Eigen::Quaternionf(
              Eigen::AngleAxisf(.5f, Eigen::Vector3f::UnitZ()))) < 1e-6f);
}

TEST_CASE("Free-joint translation uses body velocity transformed by the "
          "initial orientation") {
  using namespace rapide;
  JointCoordinates<double> joint = FreeJoint<double>{};
  Eigen::VectorXd q(7), v(6);
  const Eigen::Quaterniond orientation(
      Eigen::AngleAxisd(std::numbers::pi / 2, Vec3::UnitZ()));
  q << orientation.coeffs(), 1, 2, 3;
  v << 2, 0, 0, 1, 0, 0;
  const auto next = integration::integrate_joint_configuration(joint, q, v, .5);
  REQUIRE((next.tail<3>() - Vec3{1, 2.5, 3}).norm() < 1e-12);
  const Eigen::Quaterniond actual(next(3), next(0), next(1), next(2));
  REQUIRE(
      actual.angularDistance(orientation * Eigen::Quaterniond(Eigen::AngleAxisd(
                                               1., Vec3::UnitX()))) < 1e-12);
  JointCoordinates<double> hinge = RevoluteJoint<double>{Vec3::UnitZ()};
  q = Eigen::VectorXd::Constant(1, .3);
  v = Eigen::VectorXd::Constant(1, 2);
  REQUIRE(
      std::abs(integration::integrate_joint_configuration(hinge, q, v, .1)(0) -
               .5) < 1e-12);
}

TEST_CASE("Central Euler schemes choose old or new velocity explicitly") {
  auto advance = [](double position, double velocity, double dt) {
    return rapide::integration::advance_linear(position, velocity, dt);
  };
  const auto explicit_step =
      rapide::integration::explicit_euler(1., 2., 3., .1, advance);
  const auto semi_implicit =
      rapide::integration::semi_implicit_euler(1., 2., 3., .1, advance);
  REQUIRE(std::abs(explicit_step.configuration - 1.2) < 1e-12);
  REQUIRE(std::abs(semi_implicit.configuration - 1.23) < 1e-12);
  REQUIRE(std::abs(explicit_step.velocity - 2.3) < 1e-12);
  REQUIRE(explicit_step.velocity == semi_implicit.velocity);
}

TEST_CASE("ABA and maximal free bodies use the same finite rotation update") {
  using namespace rapide;
  ModelBuilder builder;
  const auto body = builder.add_body("body", 1, Mat3::Identity());
  auto model = builder.build();
  const Vec3 omega{1, 2, 3};
  model.initial_state.bodies[body].omega = omega;
  Simulator aba_sim(model, Backend::ABA), maximal_sim(model, Backend::Maximal);
  aba_sim.gravity.setZero();
  maximal_sim.gravity.setZero();
  aba_sim.step(.7);
  maximal_sim.step(.7);
  const auto actual = aba_sim.state().bodies[body].pose.orientation;
  REQUIRE(actual.angularDistance(
              maximal_sim.state().bodies[body].pose.orientation) < 1e-12);
  REQUIRE(actual.angularDistance(Eigen::Quaterniond(Eigen::AngleAxisd(
              .7 * omega.norm(), omega.normalized()))) < 1e-12);
}

TEST_CASE("Maximal pose projection preserves world angular velocity") {
  using namespace rapide;
  rapide::ModelBuilder builder;
  const auto body =
      builder.add_body("body", 1, Mat3::Identity(),
                       {{1, 2, 3}, rotation_increment(Vec3{.3, .4, .5})});
  const auto model = builder.build();
  auto state = model.initial_state;
  state.bodies[body].omega = Vec3{1, 2, 3};
  const Vec3 before =
      state.bodies[body].pose.orientation * state.bodies[body].omega;
  Eigen::VectorXd increment(6);
  increment << .4, -.2, .1, 1, 0, 0;
  maximal::retract(model, state, increment);
  const auto &after = state.bodies[body];
  REQUIRE((after.pose.orientation * after.omega - before).norm() < 1e-12);
  REQUIRE((after.pose.position - Vec3{2, 2, 3}).norm() < 1e-12);
}
