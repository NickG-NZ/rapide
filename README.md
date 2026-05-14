# rapide

Minimal multibody physics simulator for autonomous wheeled vehicles, written in C++.
Targets autonomous motorcycle simulation with articulated robot arms. Differentiable
by design: every algorithm is templated on `Scalar` so autodiff types (e.g.
`Eigen::AutoDiffScalar`) drop in without code changes.

## Scope

- Rigid-body dynamics via Featherstone's Articulated Body Algorithm (ABA)
- Kinematic tree representation (floating base + joints)
- Tire spring-damper contact (flat ground); no general contact solver
- No suspension or steering in v1 — just rolling chassis physics
- No external physics library dependency (Eigen only)

## Architecture

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
   assembled into a `std::vector<Force6<S>>` (one per body) before the ABA call.
   The algorithm itself is pure dynamics.

### Module map

```
include/rapide/
├── math/
│   ├── se3.hpp          SE3<S>: (R,p) rigid transform + spatial actions
│   └── spatial.hpp      Motion6<S>, Force6<S>, SpatialInertia<S>
├── model/
│   ├── joint.hpp        FreeJoint<S>, RevoluteJoint<S>; Joint<S> = std::variant
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
src/
├── algorithm/aba.cpp
├── algorithm/kinematics.cpp
└── contact/tire.cpp

vehicles/
├── motorcycle.hpp       build_motorcycle<S>() → Model<S>
└── arm.hpp              build_arm<S>(n_dof) → Model<S>

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

- **Eigen 3.4+** — all linear algebra; autodiff via `Eigen::AutoDiffScalar`
- **Catch2 v3** — unit tests
- **CMake 3.20+**

## Build

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build
```
