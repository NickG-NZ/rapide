#include <catch2/catch_test_macros.hpp>
#include "rapide/algorithm/aba.hpp"
#include "rapide/algorithm/kinematics.hpp"
#include "rapide/model/model.hpp"
#include "rapide/model/joint.hpp"
#include "rapide/model/body.hpp"
#include "vehicles/motorcycle.hpp"

using namespace rapide;
using namespace rapide::vehicles;
using S = double;
using V3 = Eigen::Matrix<S,3,1>;

static constexpr S eps = 1e-6;

// ── Single floating body free-fall ────────────────────────────────────────────
TEST_CASE("ABA: single floating body free-fall") {
    ModelBuilder<S> builder;
    Body<S> body = Body<S>::with_inertia("box", 10.0, V3::Zero(),
                                          Eigen::Matrix<S,3,3>::Identity() * 0.1);
    builder.add_body(std::move(body), FreeJoint<S>{}, 0, SE3<S>::Identity());
    Model<S> model = builder.build();

    State<S> state = model.zero_state();
    std::vector<Wrench<S>> f_ext(model.num_bodies(), Wrench<S>::Zero());
    Eigen::VectorX<S> tau = Eigen::VectorX<S>::Zero(model.nv);

    auto ddq = forward_dynamics(model, state, tau, f_ext, V3{0, 0, -9.81});

    REQUIRE(ddq.head(3).norm() < eps);              // no angular acceleration
    REQUIRE(std::abs(ddq(3)) < eps);                // no x acceleration
    REQUIRE(std::abs(ddq(4)) < eps);                // no y acceleration
    REQUIRE(std::abs(ddq(5) - (-9.81)) < eps);      // free-fall in z
}

// ── Prismatic joint: block on vertical slide ───────────────────────────────────
// A body attached to world via a vertical PrismaticJoint.
// With no spring and no tau, it should free-fall with ddq = -g.
TEST_CASE("ABA: prismatic joint free-fall") {
    ModelBuilder<S> builder;
    Body<S> body = Body<S>::with_inertia("block", 5.0, V3::Zero(),
                                          Eigen::Matrix<S,3,3>::Identity() * 0.1);
    builder.add_body(std::move(body),
                     PrismaticJoint<S>{V3{0, 0, 1}},  // slides along +Z
                     0, SE3<S>::Identity());
    Model<S> model = builder.build();

    State<S> state = model.zero_state();
    std::vector<Wrench<S>> f_ext(model.num_bodies(), Wrench<S>::Zero());
    Eigen::VectorX<S> tau = Eigen::VectorX<S>::Zero(model.nv);

    auto ddq = forward_dynamics(model, state, tau, f_ext, V3{0, 0, -9.81});

    REQUIRE(model.nq == 1);
    REQUIRE(model.nv == 1);
    // Prismatic joint along +Z, gravity -9.81 in Z → ddq = -9.81
    REQUIRE(std::abs(ddq(0) - (-9.81)) < eps);
}

// ── Prismatic joint: spring equilibrium ───────────────────────────────────────
// Block on vertical prismatic with spring: at rest, acceleration should equal
// (spring_force - gravity) / mass.
TEST_CASE("ABA: prismatic joint with spring force in tau") {
    ModelBuilder<S> builder;
    S mass = 2.0;
    Body<S> body = Body<S>::with_inertia("block", mass, V3::Zero(),
                                          Eigen::Matrix<S,3,3>::Identity() * 0.01);
    builder.add_body(std::move(body),
                     PrismaticJoint<S>{V3{0, 0, 1}},
                     0, SE3<S>::Identity());
    Model<S> model = builder.build();

    State<S> state = model.zero_state();
    std::vector<Wrench<S>> f_ext(model.num_bodies(), Wrench<S>::Zero());

    // tau = +mass*g → net force = 0 → ddq = 0
    Eigen::VectorX<S> tau(1);
    tau(0) = mass * 9.81;

    auto ddq = forward_dynamics(model, state, tau, f_ext, V3{0, 0, -9.81});
    REQUIRE(std::abs(ddq(0)) < eps);
}

// ── Motorcycle dimensionality ─────────────────────────────────────────────────
TEST_CASE("ABA: motorcycle (no suspension) dimensionality") {
    auto mm = build_motorcycle<S>();
    const auto& m = mm.model;
    // 1 FreeJoint (nv=6) + 2 RevoluteJoints (nv=1 each)
    REQUIRE(m.nv == 8);
    REQUIRE(m.nq == 9);  // 7 + 1 + 1
    REQUIRE(m.num_bodies() == 4);
}

TEST_CASE("ABA: motorcycle (no suspension) zero-state runs without NaN") {
    auto mm = build_motorcycle<S>();
    State<S> state = motorcycle_standing_state(mm.model);

    std::vector<Wrench<S>> f_ext(mm.model.num_bodies(), Wrench<S>::Zero());
    Eigen::VectorX<S> tau = Eigen::VectorX<S>::Zero(mm.model.nv);

    auto ddq = forward_dynamics(mm.model, state, tau, f_ext, V3{0,0,-9.81});
    REQUIRE(ddq.allFinite());
}

// ── Motorcycle with suspension ────────────────────────────────────────────────
TEST_CASE("ABA: motorcycle (suspended) dimensionality") {
    auto mm = build_motorcycle_suspended<S>();
    const auto& m = mm.model;
    // 1 FreeJoint (nv=6) + 2 PrismaticJoints (nv=1) + 2 RevoluteJoints (nv=1)
    REQUIRE(m.nv == 10);
    REQUIRE(m.nq == 11);  // 7 + 1 + 1 + 1 + 1
    REQUIRE(m.num_bodies() == 6);
    REQUIRE(mm.rear_fork_joint  >= 0);
    REQUIRE(mm.front_fork_joint >= 0);
}

TEST_CASE("ABA: motorcycle (suspended) zero-state runs without NaN") {
    auto mm = build_motorcycle_suspended<S>();
    State<S> state = motorcycle_standing_state(mm.model);

    std::vector<Wrench<S>> f_ext(mm.model.num_bodies(), Wrench<S>::Zero());
    Eigen::VectorX<S> tau = Eigen::VectorX<S>::Zero(mm.model.nv);

    auto ddq = forward_dynamics(mm.model, state, tau, f_ext, V3{0,0,-9.81});
    REQUIRE(ddq.allFinite());
}

TEST_CASE("Forward kinematics: motorcycle body count") {
    auto mm = build_motorcycle<S>();
    State<S> state = motorcycle_standing_state(mm.model);
    auto kd = forward_kinematics(mm.model, state);
    REQUIRE(static_cast<int>(kd.X_world.size()) == mm.model.num_bodies());
}
