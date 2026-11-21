#include "rapide/collision/coal.hpp"
#include "rapide/maximal/simulator.hpp"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using namespace rapide::maximal;
namespace collision = rapide::collision;
using Catch::Approx;

TEST_CASE("Coal adapter normalizes gaps and authored-frame witness points") {
  ModelBuilder b;
  auto ball =
      b.add_body("ball", 1, Mat3::Identity(), {{0, 0, .095}}, Vec3{.01, 0, 0});
  const auto m = b.build();
  collision::World world(m);
  world.add(0, std::make_shared<coal::Halfspace>(Vec3::UnitZ(), 0));
  world.add(ball, std::make_shared<coal::Sphere>(.1));
  auto hits = world.contacts(m.initial_state);
  REQUIRE(hits.size() == 1);
  const auto &c = hits[0];
  REQUIRE(c.gap == Approx(-.005).margin(1e-8));
  REQUIRE((c.normal - Vec3::UnitZ()).norm() < 1e-8);
  REQUIRE(world_point(m, m.initial_state, c.body_b, c.local_b).z() ==
          Approx(-.005).margin(1e-8));
  auto forces = zero_forces(m);
  collision::add_compliant_forces(m, m.initial_state, hits, forces,
                                  {1000, 0, 0, .01});
  REQUIRE(forces[ball].force.z() == Approx(5).margin(1e-8));
  REQUIRE((forces[0].force + forces[ball].force).norm() < 1e-10);
  auto state = m.initial_state;
  state.bodies[ball].pose.position.z() = .2;
  REQUIRE(world.contacts(state).empty());
}

TEST_CASE("Coal filters self and connected bodies but can include nonadjacent "
          "pairs") {
  ModelBuilder b;
  auto a = b.add_body("a", 1, Mat3::Identity());
  auto c = b.add_body("b", 1, Mat3::Identity(), {{.1, 0, 0}});
  b.add_joint(JointType::Ball, a, c, {{.1, 0, 0}});
  const auto m = b.build();
  collision::World world(m);
  world.add(a, std::make_shared<coal::Sphere>(.1));
  world.add(a, std::make_shared<coal::Sphere>(.1));
  world.add(c, std::make_shared<coal::Sphere>(.1));
  REQUIRE(world.contacts(m.initial_state).empty());
  world.include(a, c);
  REQUIRE(world.contacts(m.initial_state).size() == 2);
  world.exclude(a, c);
  REQUIRE(world.contacts(m.initial_state).empty());
}

TEST_CASE("Maximal body lands on Coal ground and can leave it") {
  ModelBuilder b;
  auto ball = b.add_body("ball", 1, Mat3::Identity(), {{0, 0, .3}});
  Simulator sim(b.build());
  collision::World world(sim.model());
  world.add(0, std::make_shared<coal::Halfspace>(Vec3::UnitZ(), 0));
  world.add(ball, std::make_shared<coal::Sphere>(.1));
  sim.set_force_callback([&](const Model &m, const State &s) {
    auto f = zero_forces(m);
    collision::add_compliant_forces(m, s, world.contacts(s), f,
                                    {10000, 100, .5, .01});
    return f;
  });
  for (int i = 0; i < 3000; ++i)
    sim.step(.0005);
  REQUIRE(sim.state().bodies[ball].pose.position.z() ==
          Approx(.1 - 9.81 / 10000).margin(1e-5));
  auto state = sim.state();
  state.bodies[ball].velocity.z() = 2;
  sim.set_state(state);
  for (int i = 0; i < 100; ++i)
    sim.step(.0005);
  REQUIRE(sim.state().bodies[ball].pose.position.z() > .15);
  REQUIRE(world.contacts(sim.state()).empty());
}
