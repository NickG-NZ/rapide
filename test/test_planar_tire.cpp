#include "examples/models/planar_wheel.hpp"
#include "rapide/simulator.hpp"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using namespace rapide;
using Catch::Approx;

TEST_CASE("Planar tire is unilateral and friction dissipates slip") {
  Tire tire;
  BodyState wheel;
  wheel.pose.position.z() = .4;
  wheel.velocity = Vec3{2, 0, -3};
  REQUIRE(tire.compute(wheel).wrench.force.norm() == 0);
  wheel.pose.position.z() = .32;
  auto response = tire.compute(wheel);
  REQUIRE(response.wrench.force.z() == Approx(420));
  REQUIRE(std::abs(response.wrench.force.x()) <=
          tire.friction * response.wrench.force.z());
  REQUIRE(response.wrench.force.dot(response.slip_velocity) <= 0);
  REQUIRE(response.wrench.torque.y() ==
          Approx(-wheel.pose.position.z() * response.wrench.force.x()));
  wheel.velocity.z() = 4;
  REQUIRE(tire.compute(wheel).wrench.force.norm() == 0);
  wheel.velocity = Vec3{2, 0, 0};
  wheel.omega.y() = 2 / wheel.pose.position.z();
  REQUIRE(tire.compute(wheel).wrench.force.x() == Approx(0).margin(1e-12));
}

TEST_CASE("Dropped wheel bounces then settles into planar rolling") {
  const auto example = rapide::examples::build_planar_wheel();
  Simulator sim(example.model);
  sim.set_force_callback([&](const Model &m, const State &s) {
    auto forces = zero_forces(m);
    forces[example.wheel] =
        example.tire.compute(s.bodies[example.wheel]).wrench;
    return forces;
  });
  bool touched = false, rebounded = false;
  for (int i = 0; i < 6000; ++i) {
    sim.step(.001);
    const auto &s = sim.state().bodies[example.wheel];
    touched |= example.tire.compute(s).compression > 0;
    rebounded |= touched && s.velocity.z() > .1;
    REQUIRE(std::abs(s.pose.position.y()) < 1e-12);
    REQUIRE(std::abs(s.velocity.y()) < 1e-12);
    REQUIRE(std::abs(s.omega.x()) + std::abs(s.omega.z()) < 1e-12);
    REQUIRE(s.pose.position.z() > .23);
  }
  REQUIRE(rebounded);
  const auto &s = sim.state().bodies[example.wheel];
  const auto response = example.tire.compute(s);
  REQUIRE(response.compression ==
          Approx(12 * 9.81 / example.tire.stiffness).margin(1e-8));
  REQUIRE(std::abs(response.slip) < 1e-9);
  REQUIRE(std::abs(s.velocity.z()) < 1e-8);
  REQUIRE(s.velocity.x() > 0);
  REQUIRE(s.velocity.x() < .65);
  REQUIRE(s.velocity.x() ==
          Approx(s.pose.position.z() * s.omega.y()).margin(1e-9));
  auto launched = sim.state();
  launched.bodies[example.wheel].velocity.z() = 2;
  sim.set_state(launched);
  bool airborne = false;
  for (int i = 0; i < 100; ++i) {
    sim.step(.001);
    const auto &wheel = sim.state().bodies[example.wheel];
    if (wheel.pose.position.z() > example.tire.radius) {
      airborne = true;
      REQUIRE(example.tire.compute(wheel).wrench.force.norm() == 0);
    }
  }
  REQUIRE(airborne);
}

TEST_CASE("Cambered tire applies contact force at the projected patch center") {
  Tire tire;
  tire.lateral_slip_stiffness = 230;
  BodyState wheel;
  wheel.pose.orientation = rotation_increment(Vec3{.4, 0, 0});
  wheel.pose.position = Vec3{.2, -.3, .29};
  wheel.velocity = Vec3{.7, -.2, -.1};
  wheel.omega = Vec3{.3, 1, -.2};
  const auto contact = tire.compute(wheel);
  REQUIRE(contact.compression > 0);
  REQUIRE(contact.contact_point.z() == Approx(0).margin(1e-12));
  REQUIRE(std::abs(contact.contact_point.y() - wheel.pose.position.y()) > .05);
  const Vec3 arm = contact.contact_point - wheel.pose.position;
  REQUIRE((contact.wrench.torque - arm.cross(contact.wrench.force)).norm() <
          1e-12);
  REQUIRE(contact.wrench.force.dot(contact.slip_velocity) <= 0);
  const double h = 1e-6;
  auto next = wheel;
  next.pose.position += h * wheel.velocity;
  next.pose.orientation =
      wheel.pose.orientation * rotation_increment(h * wheel.omega);
  const double rate =
      (tire.compute(next).compression - contact.compression) / h;
  const double expected = -tire.ground_normal.dot(
      wheel.velocity + (wheel.pose.orientation * wheel.omega).cross(arm));
  REQUIRE(rate == Approx(expected).margin(1e-6));

  const auto rotation = rotation_increment(Vec3{.3, -.2, .7});
  const Vec3 translation{.4, -.7, .2};
  auto transformed = wheel;
  transformed.pose.position = rotation * wheel.pose.position + translation;
  transformed.pose.orientation = rotation * wheel.pose.orientation;
  transformed.velocity = rotation * wheel.velocity;
  auto plane = tire;
  plane.ground_normal = rotation * tire.ground_normal;
  plane.ground_height = plane.ground_normal.dot(translation);
  const auto rotated = plane.compute(transformed);
  REQUIRE(rotated.compression == Approx(contact.compression));
  REQUIRE(
      (rotated.contact_point - rotation * contact.contact_point - translation)
          .norm() < 1e-12);
  REQUIRE((rotated.wrench.force - rotation * contact.wrench.force).norm() <
          1e-10);
  REQUIRE((rotated.wrench.torque - rotation * contact.wrench.torque).norm() <
          1e-10);
}

TEST_CASE("Horizontal axle degeneracy and loss of contact remain finite") {
  Tire tire;
  BodyState wheel;
  wheel.pose.orientation = rotation_increment(Vec3{std::acos(-1.) / 2, 0, 0});
  wheel.pose.position.z() = .05;
  const auto contact = tire.compute(wheel);
  REQUIRE(contact.compression == Approx(.02));
  REQUIRE(contact.wrench.force.allFinite());
  REQUIRE(contact.wrench.torque.norm() < 1e-10);
  wheel.pose.position.z() = .2;
  REQUIRE(tire.compute(wheel).wrench.force.norm() == 0);
}

TEST_CASE("Tire separates longitudinal and lateral slip with a shared friction "
          "circle") {
  Tire tire;
  tire.longitudinal_slip_stiffness = 100;
  tire.lateral_slip_stiffness = 200;
  BodyState wheel;
  wheel.pose.position.z() = .32;
  wheel.velocity = Vec3{.4, .1, 0};
  wheel.omega.y() = .4 / .32;
  auto response = tire.compute(wheel);
  REQUIRE(response.longitudinal_slip == Approx(0).margin(1e-12));
  REQUIRE(response.lateral_slip == Approx(.1));
  REQUIRE(response.longitudinal_force == Approx(0).margin(1e-12));
  REQUIRE(response.lateral_force == Approx(-20));
  REQUIRE(response.wrench.torque.x() == Approx(-6.4));
  wheel.velocity.y() = -.1;
  REQUIRE(tire.compute(wheel).lateral_force == Approx(20));
  tire.lateral_slip_stiffness = 0;
  REQUIRE(tire.compute(wheel).lateral_force == Approx(0).margin(1e-12));
  tire.lateral_slip_stiffness = 200;
  wheel.velocity = Vec3{10, 10, 0};
  wheel.omega.setZero();
  response = tire.compute(wheel);
  const double limit = tire.friction * response.wrench.force.z();
  REQUIRE(std::hypot(response.longitudinal_force, response.lateral_force) ==
          Approx(limit));
  REQUIRE(response.lateral_force / response.longitudinal_force == Approx(2));
  REQUIRE(response.wrench.force.dot(response.slip_velocity) < 0);
  wheel.pose.position.z() = .4;
  REQUIRE(tire.compute(wheel).lateral_force == 0);
  tire.lateral_slip_stiffness = -1;
  REQUIRE_THROWS_AS(tire.compute(wheel), std::invalid_argument);
}
