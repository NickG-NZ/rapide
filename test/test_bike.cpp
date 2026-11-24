#include "examples/models/bike.hpp"
#include "rapide/simulator.hpp"
#include <catch2/catch_test_macros.hpp>

using namespace rapide;

TEST_CASE(
    "Bike initially rolls without longitudinal slip and supports its weight") {
  const auto bike = rapide::examples::build_bike(2, .08);
  double support = 0;
  for (auto wheel : {bike.rear, bike.front}) {
    const auto response =
        bike.tire.compute(bike.model.initial_state.bodies[wheel]);
    REQUIRE(std::abs(response.longitudinal_slip) < 1e-12);
    REQUIRE(response.lateral_force < 0);
    support += response.wrench.force.z();
  }
  REQUIRE(std::abs(support - 38 * 9.81) < 1e-10);
}

TEST_CASE("Bike preserves straight rolling and responds symmetrically to "
          "lateral disturbances") {
  const auto run = [](double lateral) {
    const auto bike = rapide::examples::build_bike(2, lateral);
    Simulator sim(bike.model);
    sim.set_force_callback([&](const Model &m, const State &s) {
      return rapide::examples::bike_forces(bike, m, s);
    });
    for (int i = 0; i < 1000; ++i) {
      const auto report = sim.step(.001);
      REQUIRE(report.position_error < 1e-8);
      REQUIRE(report.orientation_error < 1e-8);
      REQUIRE(report.rank == 15);
    }
    return sim.state().bodies[bike.chassis];
  };
  const auto straight = run(0);
  REQUIRE(std::abs(straight.pose.position.y()) < 1e-10);
  REQUIRE(straight.pose.orientation.angularDistance(
              Eigen::Quaterniond::Identity()) < 1e-8);
  REQUIRE(std::abs(straight.velocity.x() - 2) < 1e-8);
  const auto left = run(.08), right = run(-.08);
  REQUIRE(std::abs(left.pose.position.x() - right.pose.position.x()) < 1e-9);
  REQUIRE(std::abs(left.pose.position.y() + right.pose.position.y()) < 1e-9);
  REQUIRE(std::abs(left.pose.position.z() - right.pose.position.z()) < 1e-9);
  REQUIRE(std::abs(left.pose.position.y()) > 1e-3);
  REQUIRE(std::abs((left.pose.orientation * Vec3::UnitY()).z()) > 1e-3);
}

TEST_CASE("Dropped bike lands, rebounds and settles without lateral motion") {
  const auto bike = rapide::examples::build_bike(2, 0, .15);
  Simulator sim(bike.model);
  for (auto wheel : {bike.rear, bike.front}) {
    REQUIRE(bike.tire.compute(sim.state().bodies[wheel]).wrench.force.norm() ==
            0);
  }
  bool landed = false, rebounded = false;
  sim.set_force_callback([&](const Model &m, const State &s) {
    return rapide::examples::bike_forces(bike, m, s);
  });
  double max_angle = 0, max_compression_change = 0;
  const double initial_length = .55;
  for (int i = 0; i < 6000; ++i) {
    const auto report = sim.step(.001);
    REQUIRE(report.position_error < 1e-8);
    REQUIRE(report.orientation_error < 1e-8);
    const double angle =
        joint_position(sim.model(), sim.state(), bike.suspension);
    max_angle = std::max(max_angle, std::abs(angle));
    const auto &frame = sim.state().bodies[bike.chassis];
    landed |= bike.tire.compute(sim.state().bodies[bike.rear]).compression > 0;
    rebounded |= landed && frame.velocity.z() > .05;
    REQUIRE(std::abs(frame.pose.position.y()) < 1e-4);
    REQUIRE(std::abs((frame.pose.orientation * Vec3::UnitY()).z()) < 1e-4);

    const double length = (world_point(sim.model(), sim.state(), bike.spring.b,
                                       bike.spring.point_b) -
                           world_point(sim.model(), sim.state(), bike.spring.a,
                                       bike.spring.point_a))
                              .norm();
    max_compression_change =
        std::max(max_compression_change, std::abs(length - initial_length));
  }
  REQUIRE(landed);
  REQUIRE(rebounded);
  REQUIRE(std::abs(sim.state().bodies[bike.chassis].velocity.z()) < 1e-3);
  REQUIRE(max_angle > .005);
  REQUIRE(max_compression_change > .001);
  REQUIRE(std::abs(joint_velocity(sim.model(), sim.state(), bike.suspension)) <
          1e-3);
  REQUIRE(std::abs(joint_position(sim.model(), sim.state(), bike.suspension)) <
          1e-3);
}
