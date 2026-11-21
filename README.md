# rapide

Multibody physics simulator for autonomous wheeled vehicles, written in C++20.
The maximal-coordinate backend supports closed linkages using analytic joint
constraints and a dense, rank-aware impulse solve. Coal supplies collision geometry
for compliant contact response. The original scalar-templated ABA backend remains
available for tree models and comparisons.

## Scope

- CPU, double-precision maximal-coordinate rigid-body dynamics
- Ball, fixed, revolute, and prismatic joints in arbitrary graphs, including loops
- Mass-weighted position/velocity projection and redundant-constraint handling
- Point spring-dampers, scalar joint springs, and torque/force actuation
- Coal broad phase, collision filtering, and compliant contact with sliding friction
- Legacy ABA motorcycle example and flat-ground tire model

Hard contact, joint limits, restitution, continuous collision detection, and
autodiff are not implemented in the maximal backend. Compliant contacts permit
penetration and need a timestep appropriate to stiffness and mass.

## Maximal-coordinate solver

Headers live under `rapide/maximal/`; the namespace is `rapide::maximal`.
Bodies are independent of joints. Add all bodies, then connect any pair with a
joint; a loop needs no tree-cutting or special loop joint.

```cpp
#include "rapide/maximal/simulator.hpp"

using namespace rapide::maximal;
ModelBuilder builder;  // Body 0 is static world.
auto link = builder.add_body("pendulum", 1.0, Mat3::Identity(), {{0.5, 0, 0}});
builder.add_joint(JointType::Revolute, 0, link, {}, {{-0.5, 0, 0}});
Simulator sim(builder.build());
sim.gravity = {0, -9.81, 0};
sim.step(0.001);
```

`add_body` accepts the authored body-origin pose and an optional local COM offset.
Runtime poses locate the **COM in world coordinates**; orientation maps body axes
to world axes. The six velocity components are `[omega_body; v_com_world]`.
Inertia is constant in body axes, about the COM. External forces and COM torques
are world-expressed. Joint frames and force attachment points use the authored
body coordinates, and revolute/prismatic axes are frame-local +Z.

`Simulator::step(dt, forces)` holds those forces constant over the step. Alternatively,
`set_force_callback` evaluates forces from the current model and state. Helpers in
`forces.hpp` apply point springs, joint efforts, and scalar joint springs. Revolute
spring coordinates are wrapped to [-pi, pi]; multi-turn springs need an unwrapped
coordinate maintained by the caller.

Each step predicts unconstrained velocity, solves bilateral constraint impulses,
integrates poses, then projects position and velocity. The solve uses mass-scaled
constraint rows and complete orthogonal decomposition, accommodating redundant rows
without silently softening joints. Initial configurations are assembled through the
same projection. `set_state` and `step` leave the simulator unchanged on failure.

`report()` exposes numerical rank, row count, projection iterations, and maximum
position, orientation, and velocity residuals. `constraint_impulses` come from the
pre-integration velocity solve; `projection_impulses` are separate numerical
corrections. Redundant constraints can have nonunique reactions even when motion
is unique. Residual tolerances are per component; angular residuals are in radians.

The orientation residual uses a principal rotation logarithm, with a chart boundary
at pi. Position projection is a local assembly method: supply a reasonable initial
configuration. Singular linkage positions and inconsistent constraints can fail
assembly. Explicit force integration and projection are first order; projection is
not energy-preserving. Dense matrices prioritize small-mechanism correctness over
large-scene throughput, and workspaces currently allocate during solves.

### Closed-loop example

`examples/models/four_bar.hpp` constructs a parallelogram with three moving bodies
and four hinges. Its 20 constraint rows have rank 17, leaving one degree of freedom.

```bash
./build/four_bar_rollout > four_bar.csv
```

The five-second driven rollout reports loop errors and rank. Tests also compare its
motion against the analytic effective inertia and check timestep convergence.

### Coal collision adapter

Link `rapide_collision` and include `rapide/collision/coal.hpp` to use the selected
collision backend. `rapide::collision::World` owns geometry and a dynamic AABB tree.
Attach a Coal primitive or convex geometry to a body with `add(body, geometry, pose)`.
Static ground can be a `coal::Halfspace` attached to body 0. Queries update shape
transforms from the supplied maximal state; they do not integrate or solve dynamics.

By default, same-body, fixed-cluster, neighboring-cluster, and static-static pairs
are filtered. `include`/`exclude` change the body-pair policy; same-body and
static-static pairs remain excluded. A `World` is mutable and cannot be shared by
concurrent queries; separate worlds can be used for independent simulations.

Contacts contain body/collider IDs, authored-local witness points, a world midpoint,
an A-to-B world normal, and a signed gap (negative penetration). Coal 3.0.3 is the
tested API. `add_compliant_forces` adds equal-and-opposite penalty forces at the
contact midpoint with regularized sliding friction. It does not provide static
sticking or a persistent multi-point manifold. Stiffness/damping are per contact.
The existing specialized tire model remains in the ABA motorcycle example.

## Architecture

The following describes the retained ABA backend. The maximal backend is separate
and uses `double`, independent body states, and graph constraints rather than
scalar templates and generalized joint coordinates.

### Design principles

1. **`Scalar` template everywhere** — `using S = double` by default; swap to
   `AutoDiffScalar<...>` or `ceres::Jet` at the call site for gradients.
2. **Store `(R, p)`, act on demand** — `SE3<S>` carries a 3×3 rotation and 3-vec
   translation. Spatial actions on 6D vectors are inline methods; the 6×6 adjoint
   matrix is never stored, only materialised via `to_adjoint()` when explicitly needed.
3. **Model ≠ State** — `Model<S>` is topology + inertia (const after build).
   `State<S>` is `(q, v)`. Algorithms are free functions over both.
4. **Joints via `std::variant`** — zero-overhead polymorphism, works cleanly with
   the scalar template. New joint types extend the variant.
5. **`f_ext` assembly before ABA** — tire forces, gravity, and arm loads are
   assembled into a `std::vector<Wrench<S>>` (one per body) before the ABA call.
   The algorithm itself is pure dynamics.

### Module map

```
rapide/
├── math/
│   ├── se3.hpp          SE3<S>: (R,p) rigid transform + spatial actions
│   └── spatial.hpp      Twist<S>, Wrench<S>, SpatialInertia<S>
├── model/
│   ├── joint.hpp        Free/Revolute/PrismaticJoint<S>, Joint<S> variant
│   ├── body.hpp         Body<S>: name + SpatialInertia
│   └── model.hpp        Model<S>: tree topology, bodies, joints, q/v offsets
├── algorithm/
│   ├── kinematics.hpp   forward_kinematics() → per-body SE3 in world frame
│   └── aba.hpp          forward_dynamics() → joint accelerations (vdot)
├── contact/
│   └── tire.hpp         TireContact<S>: spring-damper normal + linear slip
├── integrator/
│   └── integrator.hpp   rk4_step(), euler_step(); per-joint q integration
└── simulator.hpp        Simulator<S>: owns Model+State, drives time stepping
```

```
examples/
├── motorcycle_rollout.cpp
└── models/
    ├── motorcycle.hpp   build_motorcycle<S>() → Model<S>
    └── arm.hpp          build_arm<S>(n_dof) → Model<S>

test/
├── test_se3.cpp
├── test_spatial.cpp
├── test_aba.cpp
└── test_tire.cpp
```

### Kinematic tree convention

- Body index 0 is the **world** (fixed, infinite mass).
- `joints[i]` connects `bodies[i+1]` to `parent[i]`.
- `X_tree[i]` is the fixed SE3 from parent-body frame to the joint's predecessor
  frame (set at model-build time, constant).
- Joint state occupies a contiguous slice of `q` and `v`; `Model` stores per-joint
  offsets `q_idx[i]`, `v_idx[i]`.

### Floating base (FreeJoint)

`q` slice: 7 scalars `[qx, qy, qz, qw, px, py, pz]` (unit quaternion + position).  
`v` slice: 6 scalars `[ωx, ωy, ωz, vx, vy, vz]` (body-frame spatial velocity).

Configuration integration uses first-order quaternion kinematics:
```
q_new = q ⊗ exp(ω * dt / 2)   then normalise
p_new = p + R * v_body * dt
```

### v1 motorcycle topology

```
world (fixed)
└── chassis          [FreeJoint,     6 DOF]
    ├── rear_wheel   [RevoluteJoint, axis +Y]  ← drive torque
    └── front_wheel  [RevoluteJoint, axis +Y]
```

Tire contacts on rear and front wheel bodies. Gravity pre-applied as body forces.

### Planned extensions (not in v1)

| Feature             | Plug-in point                                      |
|---------------------|----------------------------------------------------|
| Suspension          | PrismaticJoint + new links                         |
| Steering fork       | RevoluteJoint (steer axis) before front wheel      |
| 3-DOF arm           | build_arm() subtree attached to chassis body       |
| Pacejka tires       | drop-in replacement for TireContact::compute()     |
| Autodiff            | instantiate Simulator<AutoDiffScalar<...>>         |
| Terrain heightmap   | extend TireContact to query height function        |

## Dependencies

### Selected collision backend for the maximal-coordinate implementation

[Coal](https://github.com/coal-library/coal) is the selected CPU collision backend
and is linked through `rapide_collision`. The retained ABA simulator still uses
its existing flat-ground tire calculation.

- Reuse Coal's broad-phase acceleration and primitive/convex geometry queries.
- Keep shape ownership, self-collision exclusions, and contact response in Rapide.
- Normalize contact normals and signed gaps at the backend boundary.
- Start with compliant tire/contact forces; hard-contact constraints are a later step.
- Add multi-point contact processing if required; do not assume a pairwise query
  supplies a persistent contact manifold.
- Target CPU execution, with multithreading available where useful; no autodiff
  requirement for this implementation.

The [collision experiment](experiments/collision/README.txt) compares Coal and
Bullet using Panda and A1 collision models, analytic contact checks, and scenes
with up to 704 shapes. See the [results and decision](experiments/collision/REPORT.txt).

### Current build dependencies

- **Eigen 3.4+** — all linear algebra; autodiff via `Eigen::AutoDiffScalar`
- **Coal 3.0.3+** — collision adapter (enabled by default)
- **Catch2 v3** — unit tests
- **CMake 3.20+**

## Build

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build
```

Use `-DCMAKE_PREFIX_PATH=/path/to/dependencies` for non-system installations.
Set `-DRAPIDE_WITH_COAL=OFF` to build the joint/loop solver without collision support.
Catch2 must be installed for CMake to register the test executables.
