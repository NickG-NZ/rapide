#include "examples/models/arm.hpp"
#include "examples/models/bike.hpp"
#include "examples/models/four_bar.hpp"
#include "rapide/simulator.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <numbers>

using namespace rapide;

TEST_CASE(
    "ABA rejects closed loops while maximal accepts the same authored model") {
  const auto linkage = examples::build_four_bar();
  REQUIRE_THROWS_WITH(
      Simulator(linkage.model, Backend::ABA),
      "ABA does not support kinematic loops; use the maximal backend");
  Simulator maximal(linkage.model, Backend::Maximal);
  REQUIRE(maximal.step(.001).rank == 17);
}

TEST_CASE("Shared free-body forces act about the COM for both backends") {
  ModelBuilder builder;
  const Pose origin{{1, 2, 3}, rotation_increment(Vec3{.2, -.4, .3})};
  const auto body =
      builder.add_body("offset", 2, Mat3::Identity(), origin, {.3, -.2, .1});
  const auto model = builder.build();
  for (Backend backend : {Backend::ABA, Backend::Maximal}) {
    Simulator sim(model, backend);
    sim.gravity.setZero();
    sim.set_force_callback([body](const Model &m, const State &) {
      auto forces = zero_forces(m);
      forces[body].force = Vec3{2, 4, -6};
      return forces;
    });
    sim.step(.01);
    const auto &s = sim.state().bodies[body];
    REQUIRE(s.omega.norm() < 1e-10);
    REQUIRE((s.velocity - Vec3{.01, .02, -.03}).norm() < 1e-10);
    REQUIRE((s.pose.position - model.initial_state.bodies[body].pose.position -
             Vec3{.0001, .0002, -.0003})
                .norm() < 1e-10);
    REQUIRE(sim.time() == .01);
  }
}

TEST_CASE(
    "Joint frame and COM conversion give the physical pendulum acceleration") {
  ModelBuilder builder;
  const auto base = builder.add_static("mount", {{.3, -.2, 1}});
  const auto body = builder.add_body(
      "pendulum", 2, Vec3{.1, .2, .3}.asDiagonal(), {{.3, -.2, 1}}, {.4, 0, 0});
  const auto hinge = rotation_increment(Vec3{-std::numbers::pi / 2, 0, 0});
  builder.add_joint(JointType::Revolute, base, body, {{0, 0, 0}, hinge},
                    {{0, 0, 0}, hinge});
  for (Backend backend : {Backend::ABA, Backend::Maximal}) {
    Simulator sim(builder.build(), backend);
    sim.step(1e-5);
    const auto &s = sim.state().bodies[body];
    REQUIRE(std::abs(s.omega.y() / 1e-5 - 2 * 9.81 * .4 / (.2 + 2 * .4 * .4)) <
            1e-6);
  }
}

TEST_CASE("ABA supports reversed graph edges, fixed joints and independent "
          "floating roots") {
  ModelBuilder builder;
  const auto child =
      builder.add_body("child", 1, Mat3::Identity(), {{1, 0, 1}});
  const auto parent =
      builder.add_body("parent", 1, Mat3::Identity(), {{0, 0, 1}});
  const auto free =
      builder.add_body("independent", 1, Mat3::Identity(), {{3, 0, 1}});
  builder.add_joint(JointType::Fixed, child, parent, {{-1, 0, 0}}, {});
  Simulator sim(builder.build(), Backend::ABA);
  sim.step(.01);
  REQUIRE((sim.state().bodies[child].pose.position -
           sim.state().bodies[parent].pose.position - Vec3{1, 0, 0})
              .norm() < 1e-10);
  REQUIRE(std::abs(sim.state().bodies[free].velocity.z() + .0981) < 1e-10);
}

TEST_CASE("Both backends accept the shared suspended bike and preserve "
          "straight rolling") {
  const auto bike = examples::build_bike(2, 0);
  Simulator aba(bike.model, Backend::ABA),
      maximal(bike.model, Backend::Maximal);
  for (auto *sim : {&aba, &maximal}) {
    sim->set_force_callback([&](const Model &m, const State &s) {
      return examples::bike_forces(bike, m, s);
    });
    for (int i = 0; i < 100; ++i) {
      sim->step(.001);
    }
  }
  for (auto body : {bike.chassis, bike.rear, bike.front, bike.swingarm}) {
    REQUIRE((aba.state().bodies[body].pose.position -
             maximal.state().bodies[body].pose.position)
                .norm() < 1e-7);
    REQUIRE((aba.state().bodies[body].velocity -
             maximal.state().bodies[body].velocity)
                .norm() < 1e-7);
  }
}

TEST_CASE("ABA rejects unsupported joints and inconsistent state assignments") {
  ModelBuilder builder;
  const auto body = builder.add_body("body", 1, Mat3::Identity());
  const auto joint = builder.add_joint(JointType::Ball, 0, body);
  REQUIRE_THROWS_AS(Simulator(builder.build(), Backend::ABA),
                    std::invalid_argument);
  auto model = builder.build();
  model.joints[joint].type = JointType::Prismatic;
  Simulator sim(model, Backend::ABA);
  auto bad = sim.state();
  bad.bodies[body].pose.position.x() = 1;
  REQUIRE_THROWS_AS(sim.set_state(bad), std::invalid_argument);
  REQUIRE(sim.state().bodies[body].pose.position.norm() == 0);
  REQUIRE(sim.time() == 0);
}

TEST_CASE("The shared simulator applies scalar efforts with "
          "backend-independent signs") {
  ModelBuilder builder;
  const auto body = builder.add_body("slider", 2, Mat3::Identity());
  const auto joint = builder.add_joint(JointType::Prismatic, 0, body);
  for (Backend backend : {Backend::ABA, Backend::Maximal}) {
    Simulator sim(builder.build(), backend);
    sim.gravity.setZero();
    sim.set_joint_effort(joint, 4);
    sim.step(.01);
    REQUIRE(std::abs(sim.state().bodies[body].velocity.z() - .02) < 1e-12);
    REQUIRE(std::abs(sim.state().bodies[body].pose.position.z() - .0002) <
            1e-12);
    const auto before = sim.state();
    REQUIRE_THROWS_AS(sim.step(-.01), std::invalid_argument);
    REQUIRE(sim.time() == .01);
    REQUIRE((sim.state().bodies[body].pose.position -
             before.bodies[body].pose.position)
                .norm() == 0);
  }
}

TEST_CASE("Two attachments to the same static world are a loop for ABA") {
  ModelBuilder builder;
  const auto mount = builder.add_static("other_mount", {{1, 0, 0}});
  const auto body = builder.add_body("bar", 1, Mat3::Identity());
  builder.add_joint(JointType::Fixed, 0, body);
  builder.add_joint(JointType::Fixed, mount, body, {{-1, 0, 0}}, {});
  REQUIRE_THROWS_WITH(
      Simulator(builder.build(), Backend::ABA),
      "ABA does not support kinematic loops; use the maximal backend");
}

TEST_CASE("The arm authoring helper produces one model for either backend") {
  ModelBuilder builder;
  const SE3<double> mount{
      rotation_increment(Vec3{.1, .2, .3}).toRotationMatrix(), {1, 2, 3}};
  examples::attach_arm<double>(builder, 0, 3, .4, 2., mount);
  const auto model = builder.build();
  REQUIRE(model.bodies.size() == 5);
  REQUIRE(model.bodies.back().mass == 5);
  for (Backend backend : {Backend::ABA, Backend::Maximal}) {
    Simulator sim(model, backend);
    sim.step(.001);
    for (const auto &body : sim.state().bodies) {
      REQUIRE(body.pose.position.allFinite());
    }
  }
}
