"""CPU collision experiment. No dynamics stepping; both backends receive identical poses.

Run with the isolated environment described in README.txt. Assets come from
pybullet_data; no visual meshes or robot dynamics are loaded into the solver.
"""
import argparse
import hashlib
import importlib.metadata
import itertools
import json
import platform
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import coal
import numpy as np
import pybullet as pb
import pybullet_data
import trimesh

RENDER = False


def transform(xyz=(0, 0, 0), rpy=(0, 0, 0)):
    t = np.eye(4)
    t[:3, :3] = np.array(pb.getMatrixFromQuaternion(pb.getQuaternionFromEuler(rpy))).reshape(3, 3)
    t[:3, 3] = xyz
    return t


def origin(element):
    e = element.find('origin')
    return transform() if e is None else transform(
        np.fromstring(e.get('xyz', '0 0 0'), sep=' '),
        np.fromstring(e.get('rpy', '0 0 0'), sep=' '))


def axis_rotation(axis, angle):
    axis = axis / np.linalg.norm(axis)
    x, y, z = axis
    k = np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]])
    t = np.eye(4)
    t[:3, :3] += np.sin(angle) * k + (1 - np.cos(angle)) * (k @ k)
    return t


class Shape:
    def __init__(self, name, link, local, geom, bullet_shape, lo, hi, radius=None, visual=-1):
        self.name, self.link, self.local = name, link, local
        self.geom = geom
        self.bullet_shape = bullet_shape
        self.obj = coal.CollisionObject(geom)
        self.body = pb.createMultiBody(baseMass=1, baseCollisionShapeIndex=bullet_shape, baseVisualShapeIndex=visual)
        pb.changeDynamics(self.body, -1, collisionMargin=0)
        self.center = (lo + hi) / 2
        self.radius = float(np.linalg.norm((hi-lo)/2) if radius is None else radius)
        self.corners = np.array(list(itertools.product(*zip(lo, hi))))

    def update(self, t):
        self.pose = t
        self.obj.setTransform(coal.Transform3s(t[:3, :3], t[:3, 3]))
        self.obj.computeAABB()
        quat = trimesh.transformations.quaternion_from_matrix(t)  # wxyz
        pb.resetBasePositionAndOrientation(self.body, t[:3, 3], np.roll(quat, -1))
        corners = self.corners @ t[:3, :3].T + t[:3, 3]
        self.lo, self.hi = corners.min(axis=0), corners.max(axis=0)
        self.world_center = t[:3, :3] @ self.center + t[:3, 3]


def make_shape(name, link, local, kind, args):
    radius = None
    visual = -1
    if kind == 'box':
        size = np.array(args)
        geom = coal.Box(*size)
        bs = pb.createCollisionShape(pb.GEOM_BOX, halfExtents=size/2)
        if RENDER:
            visual = pb.createVisualShape(pb.GEOM_BOX, halfExtents=size/2, rgbaColor=[.25,.55,.8,1])
        lo, hi = -size/2, size/2
    elif kind == 'sphere':
        r = float(args)
        geom = coal.Sphere(r)
        bs = pb.createCollisionShape(pb.GEOM_SPHERE, radius=r)
        if RENDER:
            visual = pb.createVisualShape(pb.GEOM_SPHERE, radius=r, rgbaColor=[.9,.45,.2,1])
        lo, hi, radius = np.full(3, -r), np.full(3, r), r
    elif kind == 'cylinder':
        r, length = args
        geom = coal.Cylinder(r, length)
        bs = pb.createCollisionShape(pb.GEOM_CYLINDER, radius=r, height=length)
        if RENDER:
            visual = pb.createVisualShape(pb.GEOM_CYLINDER, radius=r, length=length, rgbaColor=[.3,.7,.6,1])
        hi = np.array([r, r, length/2]); lo = -hi
        radius = np.sqrt(r*r + length*length/4)
    elif kind == 'mesh':
        path, scale = args
        mesh = trimesh.load(path, force='mesh', process=False)
        vertices = np.asarray(mesh.vertices) * scale
        # Both libraries take the convex hull of exactly the same vertex cloud.
        points = coal.StdVec_Vec3s()
        for vertex in vertices:
            points.append(vertex)
        geom = coal.Convex.convexHull(points, True, 'Qt')
        bs = pb.createCollisionShape(pb.GEOM_MESH, vertices=vertices.tolist())
        if RENDER:
            visual = pb.createVisualShape(pb.GEOM_MESH, vertices=vertices.tolist(), indices=np.asarray(mesh.faces).flatten().tolist(), rgbaColor=[.4,.6,.8,1])
        lo, hi = vertices.min(axis=0), vertices.max(axis=0)
        radius = np.linalg.norm(vertices-(lo+hi)/2, axis=1).max()
    else:
        raise ValueError(kind)
    return Shape(name, link, local, geom, bs, lo, hi, radius, visual)


def load_robot(relative):
    path = Path(pybullet_data.getDataPath()) / relative
    root = ET.parse(path).getroot()
    shapes, joints, edges, fixed_edges, hashes = [], [], set(), [], {}
    hashes[relative] = hashlib.sha256(path.read_bytes()).hexdigest()
    for link in root.findall('link'):
        for i, collision in enumerate(link.findall('collision')):
            g = list(collision.find('geometry'))[0]
            if g.tag == 'mesh':
                mesh = path.parent / g.get('filename').removeprefix('package://')
                args = mesh, np.fromstring(g.get('scale', '1 1 1'), sep=' ')
                hashes[str(mesh.relative_to(Path(pybullet_data.getDataPath())))] = hashlib.sha256(mesh.read_bytes()).hexdigest()
            elif g.tag == 'box':
                args = np.fromstring(g.get('size'), sep=' ')
            elif g.tag == 'sphere':
                args = float(g.get('radius'))
            else:
                args = float(g.get('radius')), float(g.get('length'))
            shapes.append(make_shape(f'{link.get("name")}:{i}', link.get('name'), origin(collision), g.tag, args))
    for j in root.findall('joint'):
        a, b = j.find('parent').get('link'), j.find('child').get('link')
        axis = j.find('axis')
        limits = j.find('limit')
        joints.append(dict(name=j.get('name'), kind=j.get('type'), a=a, b=b,
                           origin=origin(j), axis=np.fromstring(axis.get('xyz') if axis is not None else '1 0 0', sep=' '),
                           lower=float(limits.get('lower', '-3.14')) if limits is not None else -3.14,
                           upper=float(limits.get('upper', '3.14')) if limits is not None else 3.14))
        edges.add(frozenset((a, b)))
        if j.get('type') == 'fixed':
            fixed_edges.append((a, b))
    # Fixed-connected links are one rigid cluster. Exclude also adjacent clusters.
    group = {l.get('name'): l.get('name') for l in root.findall('link')}
    for a, b in fixed_edges:
        old, new = group[b], group[a]
        group = {k: new if v == old else v for k, v in group.items()}
    adjacent = {frozenset((group[j['a']], group[j['b']])) for j in joints}
    allowed = {(i, j) for i, a in enumerate(shapes) for j, b in enumerate(shapes) if i < j
               and group[a.link] != group[b.link]
               and frozenset((group[a.link], group[b.link])) not in adjacent}
    children = {j['b'] for j in joints}
    roots = set(group)-children
    return shapes, joints, roots, allowed, hashes


def poses(joints, roots, values, base):
    out = {r: base for r in roots}
    pending = list(joints)
    while pending:
        ready = [j for j in pending if j['a'] in out]
        if not ready:
            raise ValueError('URDF is not a tree')
        for j in ready:
            delta = np.eye(4)
            q = values.get(j['name'], 0.)
            if j['kind'] in ('revolute', 'continuous'):
                delta = axis_rotation(j['axis'], q)
            elif j['kind'] == 'prismatic':
                delta[:3, 3] = j['axis']*q
            out[j['b']] = out[j['a']] @ j['origin'] @ delta
            pending.remove(j)
    return out


def coal_contacts(a, b):
    req = coal.CollisionRequest()
    req.num_max_contacts = 8
    req.enable_contact = True
    result = coal.CollisionResult()
    coal.collide(a.obj, b.obj, req, result)
    # Coal 3.0.3 calls this penetration_depth but stores negative penetration.
    return [-float(result.getContact(i).penetration_depth) for i in range(result.numContacts())]


def sweep(shapes, allowed):
    # Conservative sweep-and-prune baseline; rebuild each frame, not an optimized tree.
    active, result = [], set()
    for i in sorted(range(len(shapes)), key=lambda i: shapes[i].lo[0]):
        a = shapes[i]
        active = [j for j in active if shapes[j].hi[0] >= a.lo[0]]
        for j in active:
            pair = tuple(sorted((i, j)))
            b = shapes[j]
            if pair in allowed and np.all(a.lo[1:] <= b.hi[1:]) and np.all(b.lo[1:] <= a.hi[1:]):
                result.add(pair)
        active.append(i)
    return result


def summary(values):
    return dict(median=float(np.median(values)), p95=float(np.percentile(values, 95)), max=float(max(values)))


def robot_case(relative, frames):
    pb.resetSimulation()
    shapes, joints, roots, allowed, hashes = load_robot(relative)
    count = len(shapes)
    ground = make_shape('ground', 'environment', np.eye(4), 'box', [8, 8, .1])
    ground.update(transform([0, 0, -.05]))
    obstacle = make_shape('obstacle', 'environment', np.eye(4), 'box', [.18, .25, .45])
    obstacle.update(transform([.35, 0, .225]))
    shapes += [ground, obstacle]
    manager = coal.DynamicAABBTreeCollisionManager()
    for shape in shapes:
        manager.registerObject(shape.obj)
    manager.setup()
    allowed |= {(i, j) for i in range(count) for j in range(count, len(shapes))}
    rng = np.random.default_rng(20261120)
    samples = [{j['name']: rng.uniform(j['lower'], j['upper']) for j in joints if j['kind'] != 'fixed'} for _ in range(16)]
    times = {k: [] for k in ('coal_all_pairs_ms', 'bullet_all_pairs_ms', 'custom_sweep_ms', 'coal_sweep_pairs_ms', 'bullet_scene_detection_ms', 'coal_scene_detection_ms')}
    candidates, sphere_counts, c_hits, b_hits = [], [], [], []
    self_aabb, self_spheres = [], []
    disagreements, misses, self_hits, env_hits = [], 0, 0, 0
    for frame in range(frames+5):
        t = max(frame-5, 0)/max(frames-1, 1)*15
        segment = min(int(t), 14); alpha = min(t-segment, 1)
        values = {k: (1-alpha)*samples[segment][k]+alpha*samples[segment+1][k] for k in samples[0]}
        base = transform([0, 0, .28 if 'a1' in relative else 0])
        world = poses(joints, roots, values, base)
        for s in shapes[:count]:
            s.update(world[s.link] @ s.local)
        start = time.perf_counter_ns(); pairs = sweep(shapes, allowed)
        sweep_ms = (time.perf_counter_ns()-start)/1e6
        spheres = {pair for pair in allowed if np.linalg.norm(shapes[pair[0]].world_center-shapes[pair[1]].world_center) <= shapes[pair[0]].radius+shapes[pair[1]].radius}
        start = time.perf_counter_ns()
        cc = {pair: coal_contacts(shapes[pair[0]], shapes[pair[1]]) for pair in sorted(allowed)}
        coal_ms = (time.perf_counter_ns()-start)/1e6
        start = time.perf_counter_ns()
        bc = {pair: pb.getClosestPoints(shapes[pair[0]].body, shapes[pair[1]].body, 0) for pair in sorted(allowed)}
        bullet_ms = (time.perf_counter_ns()-start)/1e6
        start = time.perf_counter_ns()
        for i, j in sorted(pairs):
            coal_contacts(shapes[i], shapes[j])
        filtered_ms = (time.perf_counter_ns()-start)/1e6
        start = time.perf_counter_ns(); pb.performCollisionDetection()
        scene_ms = (time.perf_counter_ns()-start)/1e6
        cb = coal.CollisionCallBackDefault()
        cb.data.request.num_max_contacts = 100000
        cb.data.request.enable_contact = True
        start = time.perf_counter_ns()
        manager.update()
        manager.collide(cb)
        coal_scene_ms = (time.perf_counter_ns()-start)/1e6
        ch = {p for p, hits in cc.items() if hits and max(hits) > 1e-5}
        bh = {p for p, hits in bc.items() if hits and min(h[8] for h in hits) < -1e-5}
        if frame < 5:  # warm-up excluded from measurements
            continue
        for key, value in zip(times, (coal_ms, bullet_ms, sweep_ms, filtered_ms, scene_ms, coal_scene_ms)):
            times[key].append(value)
        candidates.append(len(pairs)); sphere_counts.append(len(spheres))
        self_aabb.append(sum(j < count for i,j in pairs))
        self_spheres.append(sum(j < count for i,j in spheres))
        c_hits.append(len(ch)); b_hits.append(len(bh))
        misses += len((ch|bh)-pairs) + len((ch|bh)-spheres)
        self_hits += sum(j < count for i, j in ch)
        env_hits += sum(j >= count for i, j in ch)
        for pair in sorted(ch ^ bh):
            disagreements.append(dict(frame=frame, pair=[shapes[i].name for i in pair], coal_depths=cc[pair], bullet_distances=[h[8] for h in bc[pair]]))
    assert misses == 0, 'Conservative candidate filtering missed detected penetration'
    return dict(model=relative, robot_shapes=count, allowed_pairs=len(allowed), frames=frames,
                timings={k: summary(v) for k, v in times.items()}, aabb_candidates=summary(candidates),
                sphere_candidates=summary(sphere_counts), coal_colliding_pairs=summary(c_hits),
                self_aabb_candidates=summary(self_aabb), self_sphere_candidates=summary(self_spheres),
                bullet_colliding_pairs=summary(b_hits), missed_candidates=misses,
                self_collision_observations=self_hits, environment_collision_observations=env_hits,
                disagreements=disagreements, source_sha256=hashes)


def primitive_checks():
    out = []
    for kind, size, z, expected in [('sphere', .1, .095, .005), ('box', [.2, .2, .2], .095, .005), ('sphere', .1, .11, 0.)]:
        pb.resetSimulation()
        a = make_shape('object', 'object', np.eye(4), kind, size)
        b = make_shape('ground', 'ground', np.eye(4), 'box', [2, 2, .1])
        a.update(transform([0, 0, z])); b.update(transform([0, 0, -.05]))
        cc = coal_contacts(a, b)
        closest = pb.getClosestPoints(a.body, b.body, 0)
        # Persistent manifold queried separately from closest-points API.
        pb.performCollisionDetection()
        manifold = pb.getContactPoints(a.body, b.body)
        c_depth, b_depth = max(cc, default=0), max([-p[8] for p in closest], default=0)
        assert abs(c_depth-expected) < 1e-5, (kind, c_depth, expected)
        assert abs(b_depth-expected) < 1e-5, (kind, b_depth, expected)
        out.append(dict(kind=kind, height=z, expected_depth=expected, coal_depth=c_depth,
                        bullet_depth=b_depth, coal_contacts=len(cc), bullet_closest_contacts=len(closest),
                        bullet_manifold_contacts=len(manifold)))
    return out


def manifold_check():
    """Translate a slightly penetrating box smoothly along a slab without stepping."""
    pb.resetSimulation()
    a = make_shape('box', 'box', np.eye(4), 'box', [.2, .2, .2])
    b = make_shape('ground', 'ground', np.eye(4), 'box', [2, 2, .1])
    b.update(transform([0, 0, -.05]))
    counts_c, counts_b, errors_c, errors_b, normal_errors = [], [], [], [], []
    for x in np.linspace(-.1, .1, 101):
        a.update(transform([x, 0, .095]))
        req = coal.CollisionRequest(); req.enable_contact = True; req.num_max_contacts = 8
        result = coal.CollisionResult()
        coal.collide(a.obj, b.obj, req, result)
        contacts = [result.getContact(i) for i in range(result.numContacts())]
        pb.performCollisionDetection()
        bc = pb.getContactPoints(a.body, b.body)
        counts_c.append(len(contacts)); counts_b.append(len(bc))
        errors_c.extend(abs(-c.penetration_depth-.005) for c in contacts)
        errors_b.extend(abs(-c[8]-.005) for c in bc)
        # Coal's normal points A -> B; Bullet reports normal on B toward A.
        normal_errors.extend(np.linalg.norm(np.asarray(c.normal)-[0,0,-1]) for c in contacts)
        normal_errors.extend(np.linalg.norm(np.asarray(c[7])-[0,0,1]) for c in bc)
    assert max(errors_c+errors_b) < 1e-5
    assert max(normal_errors) < 1e-5
    return dict(frames=101, coal_contact_counts=summary(counts_c), bullet_contact_counts=summary(counts_b),
                max_depth_error_m=max(errors_c+errors_b), max_normal_error=max(normal_errors))


def scaling_case(copies):
    """Unfiltered native whole-scene passes; not equivalent to filtered pair queries."""
    pb.resetSimulation()
    original, joints, roots, _, _ = load_robot('a1/a1.urdf')
    values = {j['name']: (j['lower']+j['upper'])/2 for j in joints if j['kind'] != 'fixed'}
    world = poses(joints, roots, values, transform([0,0,.3]))
    shapes, initial = [], []
    for copy in range(copies):
        for source in original:
            if copy == 0:
                s = source
            else:
                s = Shape(source.name, source.link, source.local, source.geom, source.bullet_shape,
                          source.corners.min(axis=0), source.corners.max(axis=0), source.radius)
            t = world[s.link] @ s.local
            t = t.copy(); t[0,3] += 2*(copy%8); t[1,3] += 2*(copy//8)
            s.update(t); shapes.append(s); initial.append(t)
    manager = coal.DynamicAABBTreeCollisionManager()
    for s in shapes:
        manager.registerObject(s.obj)
    manager.setup()
    cm, bm = [], []
    for frame in range(35):
        for s, base in zip(shapes, initial):
            t = base.copy(); t[0,3] += .005*np.sin(frame*.1)
            s.update(t)
        cb = coal.CollisionCallBackDefault()
        cb.data.request.num_max_contacts=1000000; cb.data.request.enable_contact=True
        start=time.perf_counter_ns(); manager.update(); manager.collide(cb)
        c=(time.perf_counter_ns()-start)/1e6
        start=time.perf_counter_ns(); pb.performCollisionDetection()
        b=(time.perf_counter_ns()-start)/1e6
        if frame>=5:
            cm.append(c); bm.append(b)
    return dict(copies=copies, shapes=len(shapes), possible_pairs=len(shapes)*(len(shapes)-1)//2,
                coal_scene_ms=summary(cm), bullet_scene_ms=summary(bm), frames=30,
                final_coal_contacts=cb.data.result.numContacts(), final_bullet_contacts=len(pb.getContactPoints()))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--frames', type=int, default=300)
    parser.add_argument('--output', type=Path, default=Path(__file__).with_name('results.json'))
    args = parser.parse_args()
    assert args.frames > 0
    pb.connect(pb.DIRECT)
    result = dict(platform=platform.platform(), machine=platform.machine(),
                  packages={p: importlib.metadata.version(p) for p in ('coal', 'pybullet', 'numpy', 'trimesh')},
                  primitive_checks=primitive_checks(), manifold_check=manifold_check(), robots=[])
    for model in ('franka_panda/panda.urdf', 'a1/a1.urdf'):
        r = robot_case(model, args.frames)
        result['robots'].append(r)
        print(model, 'disagreements', len(r['disagreements']), 'misses', r['missed_candidates'], flush=True)
    result['scaling'] = [scaling_case(copies) for copies in (1, 8, 32)]
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    pb.disconnect()


if __name__ == '__main__':
    main()
