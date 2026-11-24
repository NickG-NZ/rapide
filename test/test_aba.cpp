#include "examples/models/motorcycle.hpp"
#include "rapide/aba/dynamics.hpp"
#include "rapide/aba/kinematics.hpp"
#include "rapide/aba/tree_data.hpp"
#include "rapide/model/body.hpp"
#include "rapide/model/joint.hpp"
#include <catch2/catch_test_macros.hpp>

using namespace rapide;
using namespace rapide::aba;
using namespace rapide::examples;
using S = double;
using V3 = Eigen::Matrix<S, 3, 1>;

static constexpr S eps = 1e-6;

// ── Single floating body free-fall
// ────────────────────────────────────────────
TEST_CASE("ABA: single floating body free-fall") {
  TreeBuilder<S> builder;
  SpatialBody<S> body = SpatialBody<S>::with_inertia(
      "box", 10.0, V3::Zero(), Eigen::Matrix<S, 3, 3>::Identity() * 0.1);
  builder.add_body(std::move(body), FreeJoint<S>{}, 0, SE3<S>::Identity());
  TreeData<S> model = builder.build();

  JointState<S> state = model.zero_state();
  std::vector<Wrench<S>> f_ext(model.num_bodies(), Wrench<S>::Zero());
  Eigen::VectorX<S> tau = Eigen::VectorX<S>::Zero(model.nv);

  auto ddq = forward_dynamics(model, state, tau, f_ext, V3{0, 0, -9.81});

  REQUIRE(ddq.head(3).norm() < eps);         // no angular acceleration
  REQUIRE(std::abs(ddq(3)) < eps);           // no x acceleration
  REQUIRE(std::abs(ddq(4)) < eps);           // no y acceleration
  REQUIRE(std::abs(ddq(5) - (-9.81)) < eps); // free-fall in z
}

// ── Prismatic joint: block on vertical slide
// ─────────────────────────────────── A body attached to world via a vertical
// PrismaticJoint. With no spring and no tau, it should free-fall with ddq = -g.
TEST_CASE("ABA: prismatic joint free-fall") {
  TreeBuilder<S> builder;
  SpatialBody<S> body = SpatialBody<S>::with_inertia(
      "block", 5.0, V3::Zero(), Eigen::Matrix<S, 3, 3>::Identity() * 0.1);
  builder.add_body(std::move(body),
                   PrismaticJoint<S>{V3{0, 0, 1}}, // slides along +Z
                   0, SE3<S>::Identity());
  TreeData<S> model = builder.build();

  JointState<S> state = model.zero_state();
  std::vector<Wrench<S>> f_ext(model.num_bodies(), Wrench<S>::Zero());
  Eigen::VectorX<S> tau = Eigen::VectorX<S>::Zero(model.nv);

  auto ddq = forward_dynamics(model, state, tau, f_ext, V3{0, 0, -9.81});

  REQUIRE(model.nq == 1);
  REQUIRE(model.nv == 1);
  // Prismatic joint along +Z, gravity -9.81 in Z → ddq = -9.81
  REQUIRE(std::abs(ddq(0) - (-9.81)) < eps);
}

// ── Prismatic joint: spring equilibrium
// ─────────────────────────────────────── Block on vertical prismatic with
// spring: at rest, acceleration should equal (spring_force - gravity) / mass.
TEST_CASE("ABA: prismatic joint with spring force in tau") {
  TreeBuilder<S> builder;
  S mass = 2.0;
  SpatialBody<S> body = SpatialBody<S>::with_inertia(
      "block", mass, V3::Zero(), Eigen::Matrix<S, 3, 3>::Identity() * 0.01);
  builder.add_body(std::move(body), PrismaticJoint<S>{V3{0, 0, 1}}, 0,
                   SE3<S>::Identity());
  TreeData<S> model = builder.build();

  JointState<S> state = model.zero_state();
  std::vector<Wrench<S>> f_ext(model.num_bodies(), Wrench<S>::Zero());

  // tau = +mass*g → net force = 0 → ddq = 0
  Eigen::VectorX<S> tau(1);
  tau(0) = mass * 9.81;

  auto ddq = forward_dynamics(model, state, tau, f_ext, V3{0, 0, -9.81});
  REQUIRE(std::abs(ddq(0)) < eps);
}

TEST_CASE("Shared motorcycle models run through ABA setup") {
  for (bool suspension : {false, true}) {
    MotorcycleParams<double> parameters;
    parameters.has_suspension = suspension;
    const auto bike = build_motorcycle(parameters);
    REQUIRE(bike.model.bodies.size() == (suspension ? 6 : 4));
    Simulator simulator(bike.model, Backend::ABA);
    simulator.step(.001);
    REQUIRE(simulator.state().bodies.size() == bike.model.bodies.size());
    REQUIRE(simulator.state().bodies[1].velocity.allFinite());
  }
}

TEST_CASE("ABA joint variants expose the shared joint classification") {
  REQUIRE(joint_type(JointCoordinates<double>{FreeJoint<double>{}}) ==
          JointType::Free);
  REQUIRE(joint_type(JointCoordinates<double>{
              RevoluteJoint<double>{V3::UnitY()}}) == JointType::Revolute);
  REQUIRE(joint_type(JointCoordinates<double>{
              PrismaticJoint<double>{V3::UnitZ()}}) == JointType::Prismatic);
}
