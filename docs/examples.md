# Examples

Each example is a standalone executable under `build/Release/` that writes CSV to stdout and, with `--record`, a JSON recording for the [Meshcat viewer](#meshcat-playback). Build instructions are in [development.md](development.md).

## Four-bar linkage (closed loop)

`examples/models/four_bar.hpp` constructs a parallelogram with three moving bodies
and four hinges attached to a static chassis. Its 20 constraint rows have rank 17,
leaving one degree of freedom. The builder accepts an optional mounting pose.

```bash
./build/Release/four_bar_rollout > four_bar.csv
```

The ten-second rollout mounts the mechanism in the vertical XZ plane with the
chassis fixed at z = 0.8 m. A linear spring-damper connects the chassis at local
x = 1.5 m to the coupler center (500 N/m, 50 N s/m). Its rest length is the
initial attachment separation at angle 0.9 rad. Under gravity, it compresses
and the linkage sags and settles. Meshcat draws a yellow coil and telescoping damper. The hinges enforce planar motion; velocities
are not manually clamped. This is a linkage test without tires or ground contact.

CSV output includes angle, angular velocity, coupler position, chassis height,
spring length and compression, constraint errors, and rank. With
`L = sqrt(1.25 - cos(angle))`, equilibrium satisfies
`500 * (L_rest - L) * 0.5 * sin(angle) / L = 3 * 9.81 * cos(angle)`.
Tests check that balance, fixed-chassis and planar invariants, loop closure, and
the separate torque-driven mechanism's analytic motion and timestep convergence.

## Planar wheel and deforming tire

`wheel_rollout` drops a 12 kg wheel from a hub height of 0.47 m, initially moving
at 0.65 m/s without spin. The maximal solver evolves translation and rotation;
symmetric initial conditions and forces preserve the XZ plane without velocity
clamps. This is an upright flat-ground test, not a general planar joint or a
terrain collision implementation.

`rapide/contact/compliant_tire.hpp` supplies `Tire`: an elliptic torus with outer radius
0.35 m, inner radius 0.23 m, and width 0.14 m, centered at the COM with a local-Y
axle. It supports arbitrary wheel orientation and a world-space plane
`ground_normal.dot(point) = ground_height` (unit normal required).

The lowest torus support point determines compression. Its projection onto the
plane estimates the contact-patch center. The spring-damper (4000 N/m, 100 N s/m)
uses the normal velocity at this point. Tangential slip is resolved into wheel
forward and lateral directions on the ground plane. Independent
`longitudinal_slip_stiffness` and `lateral_slip_stiffness` coefficients (both
150 N s/m by default) give `Fx = -Cx * vx_slip`, `Fy = -Cy * vy_slip` before a
shared friction-circle limit of `0.8 * normal_force`. These use slip velocities,
not slip ratio or slip angle. At an exactly vertical axle, heading is undefined;
the projected local X axis provides a deterministic fallback. There is no rolling
resistance yet: pure rolling without slip generates no tangential force. All forces act at that estimated center, with torque `arm.cross(force)`.
Slip uses the same loaded arm, so upright rolling satisfies `vx = hub_height *
omega_y`. Contact force vanishes in the air and cannot pull the tire downward.
The example uses a 1 ms timestep; stiffness and mass constrain explicit integration.

```bash
./build/Release/wheel_rollout --camber 0.35 --record /tmp/wheel.json > /tmp/wheel.csv
tools/viewer/.venv/bin/python tools/viewer/view_recording.py \
  --recording /tmp/wheel.json --html /tmp/wheel.html --open
```

Omit `--camber` for the upright drop-and-roll test. Camber is an initial angle in
radians, not a constraint: a free tilted wheel can tip and move sideways. CSV
includes the estimated patch center, compression, slip speed, signed longitudinal
and lateral slip velocities, and their corresponding contact forces.

The rigid rim and spoke markers rotate with the body. Tire vertices use the same
elliptic torus dimensions as physics, are transformed into world coordinates,
and flatten toward the ground along its normal. A quintic smootherstep joins
flat and undeformed regions with zero slope and curvature mismatch. Its clearance
band shrinks with compression; airborne geometry is unchanged. An orange marker
shows the estimated force application point while touching the ground.

This is a single-point compliant contact approximation, not a pressure-centroid
integration or rubber material model. It has no camber thrust, sidewall bulging,
or rim bottom-out forces. Deep compression can collapse parts of the visual mesh;
it is intended for modest tire deflections. The legacy ABA tire is unchanged.

The six-second recording contains 301 samples at 50 Hz. Meshcat playback stores
one tire mesh per sample and switches visibility, supporting pause and scrubbing
in standalone HTML. Tire geometry is sampled rather than interpolated; storage
grows with recording length. Dynamic vertex updates are a future viewer optimization.

## Two-wheel bike

`bike_rollout` uses a 30 kg rigid frame, two 3 kg wheels, a 2 kg rear rocker,
and a 1.3 m wheelbase. The rear wheel spins on one end of the rocker; a single
Y-axis pivot attaches the rocker to the chassis. Its opposite arm compresses
a linear spring-damper against a chassis bracket (5000 N/m, 150 N s/m).
Preload supports the initial rear load, with a 0.55 m installed spring length.
The front wheel remains on a rigid fork and steering is fixed straight ahead.
The rocker is blue, chassis orange, and spring yellow in Meshcat.
It starts at 2 m/s forward with zero lateral and vertical velocity, lifted 0.15 m
above its loaded equilibrium position. Wheel spin is matched to the loaded
rolling radius. Gravity drops it onto the ground and the suspension absorbs the landing.
Tires use 10000 N/m normal stiffness, 200 N s/m normal damping, and 300 N s/m
longitudinal and lateral slip stiffness. All six body motions remain active.

```bash
./build/Release/bike_rollout --record /tmp/bike.json > /tmp/bike.csv
tools/viewer/.venv/bin/python tools/viewer/view_recording.py \
  --recording /tmp/bike.json --html /tmp/bike.html --open
```

There is no rider, steering response, balance controller, or suspension travel stop. The
three-second drop demo stays upright by symmetry, without lateral constraints.
A safety cutoff still stops playback if lean exceeds 0.7 rad.
Frame-ground collision is not modeled. CSV logs lateral tire slips/forces,
chassis motion, lean, joint position error, suspension angle, spring length,
and spring compression. A vertical-disturbance test verifies suspension travel
and damping back to equilibrium. The initial velocities can be
changed through `build_bike(forward_speed, lateral_speed, drop_height)`.


## Meshcat playback

The optional Python viewer displays the C++ solver's recorded body poses. Meshcat
provides the browser UI and timeline; Python is not used to simulate the linkage.

```bash
python3 -m venv tools/viewer/.venv
tools/viewer/.venv/bin/python -m pip install -r tools/viewer/requirements.txt
tools/viewer/.venv/bin/python tools/viewer/view_recording.py \
  --executable ./build/Release/four_bar_rollout --open
```

Open **Animations** in Meshcat to play, pause, or scrub. Drag to orbit and scroll
to zoom. Gray is the fixed chassis, blue the crank, orange the coupler, and green
the follower. Red and white spheres mark the two independently recorded sides of
each hinge, so misalignment is visible rather than corrected by the renderer.
The animation starts paused and covers ten seconds at 50 recorded frames/second.

To keep a recording and create a standalone HTML viewer:

```bash
./build/Release/four_bar_rollout --record /tmp/four_bar.json > /tmp/four_bar.csv
tools/viewer/.venv/bin/python tools/viewer/view_recording.py \
  --recording /tmp/four_bar.json --html /tmp/four_bar.html --export-only
```

The HTML embeds its viewer, geometry, and animation; open it without a running
simulation or viewer server. Live playback serves locally while the Python command
runs; Ctrl-C stops it. `--executable` uses a temporary recording automatically.

`rapide/visualization/recording.hpp` defines the small version-1 export format:
body names/COM offsets, box visuals, joint and spring attachment points, tire
dimensions, and uniformly sampled
COM poses with xyzw quaternions. The viewer subtracts COM offsets for authored-local
geometry. It does not infer body poses from joint angles. Rendering interpolates
between recorded samples; this is playback, not live control of the simulator.

Run the viewer integration checks after building the example:

```bash
RAPIDE_FOUR_BAR_EXECUTABLE=./build/Release/four_bar_rollout \
RAPIDE_WHEEL_EXECUTABLE=./build/Release/wheel_rollout \
  tools/viewer/.venv/bin/python -m unittest discover -s tools/viewer -p 'test_*.py'
```
