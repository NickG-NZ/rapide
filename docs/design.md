# Simulator design

How the shared model, simulator, and backends fit together. For runnable demos see [examples.md](examples.md); for building see [development.md](development.md).

## Limitations

Hard contact, joint limits, restitution, continuous collision detection, and
autodiff are not implemented in the maximal backend. Compliant contacts permit
penetration and need a timestep appropriate to stiffness and mass.

## Shared simulator and model

Author a `rapide::Model` with `rapide::ModelBuilder`, then choose a backend in
`rapide::Simulator`. Both backends return the same COM body-state representation,
keep body/joint IDs stable, and accept the same world forces and COM torques.
`SimulationBackend` in `rapide/backend.hpp` defines setup, state access, state
assignment, and advancement. A custom implementation can be supplied to the
simulator as a `std::unique_ptr<SimulationBackend>`.

`Backend::Maximal` supports joint graphs with loops. `Backend::ABA` rejects loops
(including multiple attachments back to static ground), supports fixed/revolute/
prismatic joints, and automatically creates floating roots for unanchored trees.
Ball joints currently require maximal. ABA compiles its own `TreeData` and
`JointState`; these are internal numerical representations, not another authored
model. ABA setup and `set_state` require consistent joint positions/velocities;
maximal can assemble nearby states by constraint projection.

Backend setup is in `rapide/aba/backend.hpp` and `rapide/maximal/backend.hpp`.
The common interface and examples use `rapide/simulator.hpp`.
Bodies are independent of joints. Add all bodies, then connect any pair with a
joint; a loop needs no tree-cutting or special loop joint.

```cpp
#include "rapide/simulator.hpp"

using namespace rapide;
ModelBuilder builder;  // Body 0 is static world.
auto link = builder.add_body("pendulum", 1.0, Mat3::Identity(), {{0.5, 0, 0}});
builder.add_joint(JointType::Revolute, 0, link, {}, {{-0.5, 0, 0}});
Simulator sim(builder.build(), Backend::Maximal); // Or Backend::ABA for this tree.
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
`rapide/model/forces.hpp` apply point springs, joint efforts, and scalar joint springs. Revolute
spring coordinates are wrapped to [-pi, pi]; multi-turn springs need an unwrapped
coordinate maintained by the caller. `set_joint_effort(id, effort)` and
`add_joint_spring(id, stiffness, damping, rest)` add persistent scalar actuation
and springs to `step(dt)`. The explicit `step(dt, forces)` overload uses exactly
the supplied body forces. In the shared model, `nv` is the size of the packed
body-velocity vector; ABA's smaller generalized-coordinate dimensions are internal.

### Maximal solver sequencing

See the [maximal solver README](../rapide/maximal/README.md) for loop diagrams,
coordinate conventions, and the mathematical derivation.

Each maximal step predicts unconstrained velocity, solves bilateral constraint impulses,
integrates poses, then projects position and velocity. The solve uses mass-scaled
constraint rows and complete orthogonal decomposition, accommodating redundant rows
without silently softening joints. Initial configurations are assembled through the
same projection. `set_state` and `step` leave the simulator unchanged on failure.

For maximal, `report()` exposes numerical rank, row count, projection iterations, and maximum
position, orientation, and velocity residuals. `constraint_impulses` come from the
pre-integration velocity solve; `projection_impulses` are separate numerical
corrections. Redundant constraints can have nonunique reactions even when motion
is unique. Residual tolerances are per component; angular residuals are in radians.

ABA eliminates joint constraints in coordinates and does not report projection
impulses or row diagnostics.

The maximal orientation residual uses a principal rotation logarithm, with a chart boundary
at pi. Position projection is a local assembly method: supply a reasonable initial
configuration. Singular linkage positions and inconsistent constraints can fail
assembly. Explicit force integration and projection are first order; projection is
not energy-preserving. Dense matrices prioritize small-mechanism correctness over
large-scene throughput, and workspaces currently allocate during solves.

### Joints

`rapide/model/joint.hpp` owns `JointType`, the physical `Joint` definition
(body pair and attachment frames), and ABA's `JointCoordinates<S>` variant.
The authored model stores shared `Joint` objects. Its `JointConstraints` helper
in `maximal/constraints.hpp` supplies linear/angular row counts for assembly.
ABA retains typed Free/Fixed/Revolute/Prismatic coordinate operations and exposes
classification through `joint_type()`. Maximal supports Ball, Fixed, Revolute,
and Prismatic constraints; free bodies need no joint.

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

There is one authored `Model`, one body `State` interface, and one `Simulator`.
Backend-specific validation, setup, and numerical storage implement
`SimulationBackend`. Force assembly, joint effort/spring helpers, tire contact,
and visualization consume shared model and state data.

ABA retains scalar-templated spatial algebra internally. Its tree frames can differ
from authored body frames: setup rotates/translates inertia and converts initial
states, and each step converts results back to COM body state. External wrenches
are shifted from COM torque to world-origin torque before the spatial ABA call.
The public simulator uses doubles; exposing autodiff through the facade is not
implemented.

### Module map

```
rapide/
├── math/
│   ├── linalg.hpp       Shared Vec3/Mat3 aliases and scalar-generic skew matrix
│   ├── rotation.hpp     Quaternion increments/logarithms and SO(3) Jacobian
│   ├── pose.hpp         Shared quaternion Pose and validation
│   ├── se3.hpp          SE3<S>: (R,p) rigid transform + spatial actions
│   └── spatial.hpp      Twist<S>, Wrench<S>, SpatialInertia<S>
├── model/
│   ├── model.hpp        Shared Model and ModelBuilder
│   ├── state.hpp        COM body poses and velocities
│   ├── body.hpp         Mass, COM, inertia, and body-state layout
│   ├── joint.hpp        Physical Joint and coordinate representations
│   └── forces.hpp       Common forces, springs, and joint efforts
├── backend.hpp          SimulationBackend, Backend selection, options/reports
├── simulator.hpp        Shared simulator and force/effort callbacks
├── aba/
│   ├── backend.hpp      Graph validation, tree compilation, state conversion
│   ├── tree_data.hpp    Internal TreeData, JointState, and spatial inertia
│   ├── dynamics.hpp     ABA forward_dynamics()
│   ├── kinematics.hpp   Tree forward_kinematics()
│   └── configuration.hpp  Packed coordinate integration adapter
├── integration/
│   ├── configuration.hpp  Linear and quaternion updates
│   ├── joint_configuration.hpp  Joint configuration updates
│   └── euler.hpp        Shared Euler schemes
├── maximal/
│   ├── backend.hpp      Bilateral constraint solver
│   ├── constraints.hpp  JointConstraints and Jacobian assembly
│   └── state.hpp        Body velocity packing and retraction
├── contact/
│   ├── compliant_tire.hpp  World-space tire contact for either backend
│   └── tire.hpp         Spatial-wrench tire model used by motorcycle example
└── visualization/
    └── recording.hpp    Shared body-state recording export
```


```
examples/
├── bike_rollout.cpp       Two-wheel bike with rocker suspension
├── four_bar_rollout.cpp   Closed-loop parallelogram linkage
├── motorcycle_rollout.cpp ABA motorcycle with spatial-wrench tires
├── wheel_rollout.cpp      Single deforming tire drop-and-roll
└── models/
    ├── arm.hpp
    ├── bike.hpp
    ├── four_bar.hpp
    ├── motorcycle.hpp
    └── planar_wheel.hpp

test/                      Catch2 unit tests (one file per module)
tools/viewer/              Meshcat playback of recorded rollouts
experiments/collision/     Coal vs. Bullet collision benchmark
```

### Internal ABA tree convention

- Body index 0 is the **world** (fixed, infinite mass).
- `joints[i]` connects `bodies[i+1]` to `parent[i]`.
- `X_tree[i]` is the fixed SE3 from parent-body frame to the joint's predecessor
  frame (set at model-build time, constant).
- Joint state occupies a contiguous slice of `q` and `v`; `TreeData` stores per-joint
  offsets `q_idx[i]`, `v_idx[i]`.

### Floating base (FreeJoint)

`q` slice: 7 scalars `[qx, qy, qz, qw, px, py, pz]` (unit quaternion + position).  
`v` slice: 6 scalars `[ωx, ωy, ωz, vx, vy, vz]` (body-frame spatial velocity).

Configuration updates live under `rapide/integration`. Joint types describe
kinematics but do not implement integration. ABA maps packed state through the
shared joint update and semi-implicit Euler scheme. Maximal uses the same linear
updates and quaternion rotation composition, retaining velocity impulses and
pose/velocity constraint projection in its solver.

Orientation integration uses a body-axis exponential delta quaternion (exact for
angular velocity held constant over a step). Position integration remains Euler;
this does not make the entire dynamics step exact:
```
q_new = normalize(q ⊗ rotation_increment(ω_body * dt))
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


## Collision backend selection

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

The [collision experiment](../experiments/collision/README.txt) compares Coal and
Bullet using Panda and A1 collision models, analytic contact checks, and scenes
with up to 704 shapes. See the [results and decision](../experiments/collision/REPORT.txt).

## Performance work to follow

- Benchmark CPU step latency (median and tail), allocations per step, and time
  spent assembling constraints, factorizing, projecting, and querying collisions.
  Use the four-bar and larger linkages, with and without contact, after warm-up.
- Preallocate per-simulator workspaces for residuals, Jacobians, inverse-mass
  responses, constraint matrices, factorization scratch, impulses, and trial
  states. Body count and bilateral row count are known when the model is built;
  reuse that capacity across steps and projection/line-search iterations.
- Verify allocation counts, including Eigen temporaries and decomposition internals,
  rather than assuming that pre-sizing the output matrices eliminates allocation.
- Reuse storage and constant mass/topology data, not numerical factorizations:
  Jacobian values and numerical rank can change with configuration. Contact counts
  also vary, so reserve capacity with explicit growth handling for contact buffers.
- Compare performance and physics residuals before and after optimization. Establish
  the allocation baseline before choosing sparse matrices or parallel execution.
