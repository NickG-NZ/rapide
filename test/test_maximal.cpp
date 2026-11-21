#include "examples/models/four_bar.hpp"
#include "rapide/maximal/simulator.hpp"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <limits>

using namespace rapide::maximal;
using Catch::Approx;

TEST_CASE("Maximal free body follows gravity and body-axis torque") {
  ModelBuilder b;
  const auto orientation = rotation_increment(Vec3{.3, .7, -.2});
  const auto id = b.add_body("body", 2, Vec3{1, 2, 3}.asDiagonal(),
                             {{1, 2, 3}, orientation}, Vec3{.1, .2, 0});
  Simulator sim(b.build());
  const Vec3 initial = sim.state().bodies[id].pose.position;
  auto f = zero_forces(sim.model());
  f[id].torque = orientation * Vec3{0, 2, 0};
  sim.step(.01, f);
  REQUIRE((sim.state().bodies[id].omega - Vec3{0, .01, 0}).norm() < 1e-12);
  REQUIRE((sim.state().bodies[id].velocity - Vec3{0, 0, -.0981}).norm() <
          1e-12);
  REQUIRE(
      (sim.state().bodies[id].pose.position - initial - Vec3{0, 0, -.000981})
          .norm() < 1e-12);
  REQUIRE(sim.state().bodies[id].pose.orientation.norm() == Approx(1));
  REQUIRE(sim.report().constraint_rows == 0);
}

TEST_CASE("Maximal joint Jacobians match tangent finite differences") {
  for (auto type : {JointType::Ball, JointType::Fixed, JointType::Revolute,
                    JointType::Prismatic}) {
    INFO("Joint type " << static_cast<int>(type));
    ModelBuilder b;
    auto a = b.add_body("a", 2, Mat3::Identity(),
                        {{.1, -.2, .3}, rotation_increment(Vec3{.2, .3, -.1})},
                        Vec3{.1, 0, -.1});
    auto c = b.add_body("b", 3, Mat3::Identity(),
                        {{.4, .3, -.2}, rotation_increment(Vec3{-.3, .1, .4})},
                        Vec3{0, .2, 0});
    b.add_joint(type, a, c,
                {{.2, -.1, .3}, rotation_increment(Vec3{.1, -.2, .1})},
                {{-.1, .2, .1}, rotation_increment(Vec3{.2, .1, -.1})});
    const auto m = b.build();
    const auto constraints = assemble_constraints(m, m.initial_state);
    Eigen::MatrixXd numeric(constraints.residual.size(), m.nv);
    const double h = 1e-6;
    for (int k = 0; k < m.nv; ++k) {
      Eigen::VectorXd increment = Eigen::VectorXd::Zero(m.nv);
      increment(k) = h;
      auto plus = m.initial_state, minus = m.initial_state;
      retract(m, plus, increment);
      retract(m, minus, -increment);
      numeric.col(k) = (assemble_constraints(m, plus).residual -
                        assemble_constraints(m, minus).residual) /
                       (2 * h);
    }
    REQUIRE((numeric - constraints.jacobian).cwiseAbs().maxCoeff() < 2e-9);
  }
}

TEST_CASE("Maximal free rotation includes the gyroscopic torque") {
  ModelBuilder b;
  auto id = b.add_body("rotor", 1, Vec3{1, 2, 3}.asDiagonal());
  Simulator sim(b.build());
  sim.gravity.setZero();
  auto state = sim.state();
  state.bodies[id].omega = Vec3{1, 2, 3};
  sim.set_state(state);
  sim.step(.001);
  REQUIRE(
      (sim.state().bodies[id].omega - Vec3{.994, 2.003, 3 - .002 / 3}).norm() <
      1e-12);
}

TEST_CASE("Maximal fixed joint assembles and supports an offset COM") {
  ModelBuilder b;
  const Vec3 com{.3, 0, 0};
  auto id = b.add_body(
      "fixed", 2, Mat3::Identity(),
      {{.01, -.02, .01}, rotation_increment(Vec3{.01, .02, -.03})}, com);
  b.add_joint(JointType::Fixed, 0, id);
  Simulator sim(b.build());
  for (int i = 0; i < 100; ++i)
    sim.step(.01);
  REQUIRE((sim.state().bodies[id].pose.position - com).norm() < 1e-9);
  REQUIRE(sim.state().bodies[id].velocity.norm() < 1e-10);
  REQUIRE(sim.report().rank == 6);
  REQUIRE(sim.report().constraint_impulses(2) == Approx(.1962).margin(1e-9));
}

TEST_CASE("Maximal slider supports spring equilibrium and actuation") {
  ModelBuilder b;
  auto id = b.add_body("slider", 2, Mat3::Identity(), {{0, 0, .4}});
  auto j = b.add_joint(JointType::Prismatic, 0, id);
  Simulator sim(b.build());
  sim.set_force_callback([j](const Model &m, const State &s) {
    auto f = zero_forces(m);
    add_joint_spring(m, s, f, j, 100, 5, .5962);
    return f;
  });
  for (int i = 0; i < 100; ++i)
    sim.step(.002);
  REQUIRE(joint_position(sim.model(), sim.state(), j) ==
          Approx(.4).margin(1e-10));
  REQUIRE(sim.report().rank == 5);
  auto f = zero_forces(sim.model());
  add_joint_effort(sim.model(), sim.state(), f, j, 21.62);
  sim.step(.01, f);
  REQUIRE(joint_velocity(sim.model(), sim.state(), j) ==
          Approx(.01).margin(1e-10));
}

TEST_CASE("Maximal ball joint permits rotation while holding its anchor") {
  ModelBuilder b;
  auto id = b.add_body("pendulum", 1, Mat3::Identity(), {{1, 0, 0}});
  b.add_joint(JointType::Ball, 0, id, {}, {{-1, 0, 0}});
  Simulator sim(b.build());
  for (int i = 0; i < 1000; ++i)
    sim.step(.001);
  REQUIRE(world_point(sim.model(), sim.state(), id, Vec3{-1, 0, 0}).norm() <
          2e-9);
  REQUIRE(sim.state().bodies[id].pose.position.z() < -.1);
  REQUIRE(sim.report().rank == 3);
}

TEST_CASE("Maximal four-bar handles redundant rows and retains one degree of "
          "freedom") {
  const auto mechanism = rapide::examples::build_four_bar();
  Simulator sim(mechanism.model);
  sim.gravity.setZero();
  REQUIRE(sim.report().constraint_rows == 20);
  REQUIRE(sim.report().rank == 17);
  auto state = sim.state();
  state.bodies[mechanism.coupler].pose.position += Vec3{.005, -.003, .002};
  sim.set_state(state);
  const double initial =
      joint_position(sim.model(), sim.state(), mechanism.drive);
  for (int i = 0; i < 2000; ++i) {
    auto f = zero_forces(sim.model());
    add_joint_effort(sim.model(), sim.state(), f, mechanism.drive, .3);
    const auto report = sim.step(.001, f);
    REQUIRE(report.position_error < 1e-9);
    REQUIRE(report.orientation_error < 1e-9);
    REQUIRE(report.velocity_error < 1e-8);
  }
  REQUIRE(joint_position(sim.model(), sim.state(), mechanism.drive) >
          initial + .05);
  REQUIRE(sim.report().rank == 17);
}

TEST_CASE("Maximal spring forces conserve resultant force and moment") {
  ModelBuilder b;
  auto a = b.add_body("a", 1, Mat3::Identity(), {{0, 0, 0}});
  auto c = b.add_body("b", 1, Mat3::Identity(), {{1, 1, 0}});
  const auto m = b.build();
  auto f = zero_forces(m);
  add_spring(m, m.initial_state, f, a, Vec3{0, .2, 0}, c, Vec3{.1, 0, 0}, 100,
             2, .5);
  REQUIRE((f[a].force + f[c].force).norm() < 1e-12);
  const Vec3 moment = f[a].torque + f[c].torque +
                      m.initial_state.bodies[c].pose.position.cross(f[c].force);
  REQUIRE(moment.norm() < 1e-12);
}

TEST_CASE("Maximal four-bar converges to its analytic driven motion") {
  const auto mechanism = rapide::examples::build_four_bar();
  const auto rollout = [&](double dt) {
    Simulator sim(mechanism.model);
    sim.gravity.setZero();
    for (int i = 0; i < static_cast<int>(std::round(.5 / dt)); ++i) {
      auto f = zero_forces(sim.model());
      add_joint_effort(sim.model(), sim.state(), f, mechanism.drive, .3);
      sim.step(dt, f);
    }
    return joint_position(sim.model(), sim.state(), mechanism.drive);
  };
  const double effective_inertia = 2 * (1. / 12 + .25) + 2;
  const double expected = .6 + .5 * (.3 / effective_inertia) * .5 * .5;
  const double coarse = std::abs(rollout(.002) - expected);
  const double fine = std::abs(rollout(.0005) - expected);
  REQUIRE(fine < coarse * .4);
  REQUIRE(fine < 2e-5);
}

TEST_CASE("Maximal hinge rotates through full turns without locking") {
  ModelBuilder b;
  const auto id = b.add_body("rotor", 1, Mat3::Identity());
  b.add_joint(JointType::Revolute, 0, id);
  Simulator sim(b.build());
  sim.gravity.setZero();
  auto state = sim.state();
  state.bodies[id].omega.z() = 4;
  sim.set_state(state);
  for (int i = 0; i < 2000; ++i)
    sim.step(.002);
  REQUIRE(sim.state().bodies[id].omega.z() == Approx(4).margin(1e-9));
  REQUIRE(sim.state().bodies[id].pose.orientation.angularDistance(
              rotation_increment(Vec3{0, 0, 16})) < 1e-9);
  REQUIRE(sim.report().rank == 5);
}

TEST_CASE("Maximal failure leaves state and time unchanged") {
  const auto mechanism = rapide::examples::build_four_bar();
  Simulator sim(mechanism.model);
  const auto before = sim.state();
  REQUIRE_THROWS_AS(sim.step(-.01), std::invalid_argument);
  auto forces = zero_forces(sim.model());
  forces[mechanism.crank].force.x() = std::numeric_limits<double>::quiet_NaN();
  REQUIRE_THROWS_AS(sim.step(.01, forces), std::invalid_argument);
  REQUIRE(sim.time() == 0);
  REQUIRE((sim.state().bodies[mechanism.crank].pose.position -
           before.bodies[mechanism.crank].pose.position)
              .norm() == 0);
  auto invalid = before;
  invalid.bodies[0].pose.position.x() = 1;
  REQUIRE_THROWS_AS(sim.set_state(invalid), std::invalid_argument);
}

TEST_CASE("Maximal inconsistent constraints report assembly failure") {
  ModelBuilder b;
  auto id = b.add_body("body", 1, Mat3::Identity());
  b.add_joint(JointType::Ball, 0, id);
  b.add_joint(JointType::Ball, 0, id, {{1, 0, 0}});
  REQUIRE_THROWS_AS(Simulator(b.build()), std::runtime_error);
  REQUIRE_THROWS_AS(b.add_body("invalid", 1, Mat3::Zero()),
                    std::invalid_argument);
}
