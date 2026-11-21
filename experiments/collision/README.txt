Rapide CPU collision-library experiment
======================================

Purpose
-------
Evaluate Coal and Bullet as geometry/contact providers for Rapide's planned
maximal-coordinate solver. This experiment does not modify Rapide's solver.
It does not impose a single-thread restriction on the eventual implementation.

Run
---
Create a separate Python environment, install requirements.txt, then run:

    python experiments/collision/compare.py --frames 300
    python experiments/collision/render_models.py

Results are written to experiments/collision/results.json. Timings include the
Python API costs indicated below; they are not a native C++ throughput ranking.
The run prints the number of disagreements and missed broad-phase candidates.

For the initial local experiment, the environment is:
    /tmp/rapide-collision-venv/bin/python

On this Mac, PyBullet 3.2.7 needed a source build and a platform compatibility
patch in examples/ThirdPartyLibs/zlib/zutil.h. The condition
    #if defined(MACOS) || defined(TARGET_OS_MAC)
was replaced with
    #if (defined(MACOS) || defined(TARGET_OS_MAC)) && !defined(__APPLE__)
to avoid redefining fdopen for modern macOS. This is the bundled compression
library's classic-Mac compatibility branch, not collision code. The patched
source is in /tmp/rapide-collision-downloads/pybullet-3.2.7. Build with pip install
on that directory. The initial unsuccessful build logs are also under /tmp.

Models and provenance
---------------------
Both models are loaded from the pybullet_data package bundled with PyBullet
3.2.7; source hashes for URDFs and collision meshes are recorded in results.json.

  Franka Panda: franka_panda/panda.urdf
  Unitree A1:   a1/a1.urdf

Upstream:
https://github.com/bulletphysics/bullet3/tree/master/examples/pybullet/gym/pybullet_data
https://github.com/bulletphysics/bullet3/tree/master/examples/pybullet/gym/pybullet_data/franka_panda
https://github.com/bulletphysics/bullet3/tree/master/examples/pybullet/gym/pybullet_data/a1

Assets are not copied into this repository. Refer to each upstream asset's
license before redistributing it; library licensing does not automatically
cover all bundled model assets.

Method
------
- Parse URDF collision geometry, origins, joints, and joint limits.
- Compute forward kinematics once in Python and give identical transforms to
  both collision libraries. Bullet does NOT integrate a multibody or solve joints.
- Use identical primitive dimensions and identical input vertex clouds for
  convex mesh hulls. This intentionally convexifies each Panda collision mesh;
  it is not a concave triangle-mesh comparison.
- Set Bullet collision margins to zero; treat penetration above 10 micrometres
  as collision for disagreement counts. Raw depths are retained on disagreements.
- Exclude same-rigid-cluster and directly connected rigid-cluster pairs for
  explicit self-collision queries. This is a simple policy, not a robot-specific
  allowed-collision matrix. Ground/obstacle collisions remain eligible.
- Sweep 300 configurations interpolated between 16 seeded joint-limit samples.
  These are stress configurations, not validated robot motions or trajectories.
- Use an 8 x 8 m ground slab and a fixed box obstacle. No physical stepping.
- Compare exhaustive Coal contacts and Bullet getClosestPoints for the same
  eligible pairs. Independently check conservative world-AABB sweep-and-prune
  and single bounding spheres against those detected penetrations.
- Measure Coal narrow phase after custom AABB filtering separately from the
  Python filtering cost.
- Also measure unfiltered native whole-scene passes: Coal's dynamic AABB tree
  plus default collision callback, and Bullet performCollisionDetection with
  its persistent contact cache. These have different contact-generation policies
  and include connected-body collisions. They are NOT the same task as the
  filtered exhaustive pair benchmark. Coal timing includes manager.update;
  common transform updates and Python result extraction are excluded.
- Scale native passes to 1, 8, and 32 spatially separated A1 instances, translated
  coherently. This tests sparse scenes and temporal coherence, not dense piles.
- Validate sphere/slab and box/slab penetration against analytic answers; check
  normal direction, depth, and contact count while sliding a box through 101 poses.

Limitations
-----------
No hard-contact solve, dynamics stability test, continuous collision detection,
high-speed tunneling test, concave terrain, capsule fitting, or C++ benchmark.
The candidate-count test evaluates a single sphere per shape, not an optimized
multi-sphere or capsule approximation. The custom broad phase is a Python
sweep-and-prune baseline, not a persistent hierarchical BVH implementation.
No real-time capability or collision completeness beyond these samples is claimed.
The test should guide a backend prototype, not establish universal speed rankings.
models.png is a visual check of the source collision geometry: Panda on the left,
A1 on the right. Panda rendering uses the source mesh triangles; queries use their
convex hulls. It is not a screenshot of an integrated Rapide simulator.

Backend references
------------------
Coal:   https://github.com/coal-library/coal
Bullet: https://github.com/bulletphysics/bullet3
