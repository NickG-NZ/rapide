"""Display recorded Rapide body poses in Meshcat; simulation stays in C++."""
import argparse
import json
from pathlib import Path
import subprocess
import tempfile
import time

import meshcat
import meshcat.geometry as geometry
from meshcat.animation import Animation
from meshcat import transformations
import numpy as np
from tire_geometry import populate_tires


def matrix(pose):
    """Convert world COM position and xyzw quaternion into a rigid transform."""
    pose = np.asarray(pose, dtype=float)
    if pose.shape != (7,) or not np.isfinite(pose).all():
        raise ValueError("Invalid recorded pose")
    if abs(np.linalg.norm(pose[3:]) - 1) > 1e-7:
        raise ValueError("Recorded quaternion is not normalized")
    result = transformations.quaternion_matrix(pose[[6, 3, 4, 5]])
    result[:3, 3] = pose[:3]
    return result


def load_recording(path):
    """Validate the version and uniform sample times used by Meshcat animation."""
    scene = json.loads(path.read_text())
    if scene.get("version") != 1:
        raise ValueError("Unsupported recording version")
    fps = scene["fps"]
    if not np.isfinite(fps) or fps <= 0 or not scene["frames"]:
        raise ValueError("Recording needs a positive rate and at least one frame")
    for index, frame in enumerate(scene["frames"]):
        if not np.isfinite(frame["time"]) or abs(frame["time"] - index / fps) > 1e-6:
            raise ValueError("Frames must be uniformly sampled, starting at time zero")
        if len(frame["poses"]) != len(scene["bodies"]):
            raise ValueError("Frame body count does not match scene")
        for pose in frame["poses"]:
            matrix(pose)
    for body in scene["bodies"]:
        com = np.asarray(body["com"], dtype=float)
        if com.shape != (3,) or not np.isfinite(com).all():
            raise ValueError("Invalid COM offset")
    return scene


def spring_pose(scene, frame, spring):
    """Return a rigid Y-axis frame and length between the physical attachment points."""
    points = []
    for side in ("a", "b"):
        body = spring[side]
        local = np.asarray(spring[f"point_{side}"]) - scene["bodies"][body]["com"]
        points.append((matrix(frame["poses"][body]) @ np.r_[local, 1])[:3])
    direction = points[1] - points[0]
    length = np.linalg.norm(direction)
    if length < 1e-9:
        raise ValueError("Cannot draw a zero-length spring")
    y = direction / length
    reference = np.eye(3)[np.argmin(np.abs(y))]
    x = np.cross(y, reference)
    x /= np.linalg.norm(x)
    transform = np.eye(4)
    transform[:3, :3] = np.column_stack((x, y, np.cross(x, y)))
    transform[:3, 3] = points[0]
    return transform, length


def coil_geometry():
    """Build a tubular eight-turn helix with normalized axial length."""
    phase = np.linspace(0, 16 * np.pi, 257)
    centers = np.column_stack((.065 * np.cos(phase), np.linspace(.12, .88, len(phase)), .065 * np.sin(phase)))
    tangent = np.column_stack((-.065 * np.sin(phase), np.full(len(phase), .76 / (16 * np.pi)), .065 * np.cos(phase)))
    tangent /= np.linalg.norm(tangent, axis=1)[:, None]
    radial = np.column_stack((np.cos(phase), np.zeros(len(phase)), np.sin(phase)))
    normal = np.cross(tangent, radial)
    ring = np.arange(8) * 2 * np.pi / 8
    vertices = centers[:, None, :] + .009 * (np.cos(ring)[None, :, None] * radial[:, None, :] + np.sin(ring)[None, :, None] * normal[:, None, :])
    faces = []
    for i in range(len(phase) - 1):
        for j in range(8):
            a, b = i * 8 + j, i * 8 + (j + 1) % 8
            faces.extend(((a, b, a + 8), (b, b + 8, a + 8)))
    return geometry.TriangularMeshGeometry(vertices.reshape(-1, 3), faces)


def frame_camera(viewer, scene):
    """Frame moving bodies across playback, excluding world and ground visuals."""
    poses = np.asarray([frame["poses"] for frame in scene["frames"]])
    points = poses[:, 1:, :3].reshape(-1, 3)
    if not len(points):
        return
    center = (points.min(axis=0) + points.max(axis=0)) / 2
    distance = max(3., np.linalg.norm(np.ptp(points, axis=0)) + 2.)
    eye = distance * np.array([.35, .3, 1.])
    back = eye / np.linalg.norm(eye)
    right = np.cross([0., 1., 0.], back)
    right /= np.linalg.norm(right)
    transform = np.eye(4)
    transform[:3, :3] = np.column_stack((right, np.cross(back, right), back))
    quat = transformations.quaternion_from_matrix(transform)
    viewer["/Cameras/default"].set_transform(transformations.translation_matrix(center))
    viewer["/Cameras/default/rotated"].set_transform(transformations.rotation_matrix(np.pi / 2, [1, 0, 0]))
    camera = viewer["/Cameras/default/rotated/<object>"]
    camera.set_property("position", eye.tolist())
    camera.set_property("quaternion", quat[[1, 2, 3, 0]].tolist())
    camera.set_property("zoom", 1.)


def populate(viewer, scene):
    """Attach authored geometry to COM nodes and animate only recorded poses."""
    root = viewer["rapide"]
    root.delete()
    viewer["/meshcat"].set_property("visible", True)
    root.set_property("visible", True)
    frame_camera(viewer, scene)
    nodes = [root[f"body_{i}_{body['name'].replace('/', '_')}"] for i, body in enumerate(scene["bodies"])]
    for index, visual in enumerate(scene["visuals"]):
        if visual["type"] != "box":
            raise ValueError("Unsupported visual shape")
        body = visual["body"]
        local = matrix(visual["pose"])
        local[:3, 3] -= scene["bodies"][body]["com"]
        node = nodes[body][f"shape_{index}"]
        node.set_object(geometry.Box(visual["size"]), geometry.MeshLambertMaterial(color=visual["color"]))
        node.set_transform(local)
    for index, joint in enumerate(scene["joints"]):
        for side, color, radius in (("a", 0xDD4455, .045), ("b", 0xF5F5F5, .032)):
            body = joint[side]
            point = np.asarray(joint[f"point_{side}"]) - scene["bodies"][body]["com"]
            node = nodes[body][f"joint_{index}_{side}"]
            node.set_object(geometry.Sphere(radius), geometry.MeshLambertMaterial(color=color))
            node.set_transform(transformations.translation_matrix(point))
    for node, pose in zip(nodes, scene["frames"][0]["poses"]):
        node.set_transform(matrix(pose))
    shocks = []
    for index, spring in enumerate(scene.get("springs", [])):
        node = root[f"spring_{index}"]
        coil = node["coil"]
        coil.set_object(coil_geometry(), geometry.MeshLambertMaterial(color=0xFFCC22))
        housing = node["housing"]
        housing.set_object(geometry.Cylinder(.38, .028), geometry.MeshLambertMaterial(color=0x505866))
        housing.set_transform(transformations.translation_matrix([0, .19, 0]))
        rod = node["rod"]
        rod.set_object(geometry.Cylinder(.48, .013), geometry.MeshLambertMaterial(color=0xDDE5EE))
        for side in ("a", "b"):
            body = spring[side]
            point = np.asarray(spring[f"point_{side}"]) - scene["bodies"][body]["com"]
            anchor = nodes[body][f"spring_{index}_{side}"]
            anchor.set_object(geometry.Sphere(.045), geometry.MeshLambertMaterial(color=0xFFCC22))
            anchor.set_transform(transformations.translation_matrix(point))
        transform, length = spring_pose(scene, scene["frames"][0], spring)
        node.set_transform(transform)
        coil.set_property("scale", [1, length, 1])
        rod.set_transform(transformations.translation_matrix([0, length - .24, 0]))
        shocks.append((spring, node, coil, rod))
    animation = Animation(default_framerate=scene["fps"])
    for index, frame in enumerate(scene["frames"]):
        for node, pose in zip(nodes, frame["poses"]):
            with animation.at_frame(node, index) as animated:
                animated.set_transform(matrix(pose))
        for spring, node, coil, rod in shocks:
            transform, length = spring_pose(scene, frame, spring)
            with animation.at_frame(node, index) as animated:
                animated.set_transform(transform)
            with animation.at_frame(coil, index) as animated:
                animated.set_property("scale", "vector3", [1, length, 1])
            with animation.at_frame(rod, index) as animated:
                animated.set_transform(transformations.translation_matrix([0, length - .24, 0]))
    populate_tires(root, nodes, scene, animation, matrix)
    viewer.set_animation(animation, play=False, repetitions=1)
    return animation


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--recording", type=Path, help="Load an existing Rapide JSON recording")
    source.add_argument("--executable", type=Path, help="Run a rollout executable to produce a recording")
    parser.add_argument("--open", action="store_true", help="Open the viewer in the default browser")
    parser.add_argument("--html", type=Path, help="Also save a standalone interactive HTML recording")
    parser.add_argument("--export-only", action="store_true", help="Exit after saving --html")
    args = parser.parse_args()
    if args.export_only and args.html is None:
        parser.error("--export-only requires --html")
    if args.executable:
        with tempfile.TemporaryDirectory(prefix="rapide-viewer-") as directory:
            path = Path(directory) / "recording.json"
            subprocess.run([str(args.executable.resolve()), "--record", str(path)], stdout=subprocess.DEVNULL, check=True)
            scene = load_recording(path)
    else:
        scene = load_recording(args.recording)
    viewer = meshcat.Visualizer()
    populate(viewer, scene)
    if args.html:
        args.html.write_text(viewer.static_html())
        print(f"Saved {args.html.resolve()}", flush=True)
    print(f"Viewer: {viewer.url()}", flush=True)
    print(f"Loaded {len(scene['frames'])} frames, {scene['frames'][-1]['time']:.3f} seconds. Open Animations to play or scrub.", flush=True)
    if args.open:
        viewer.open()
    if not args.export_only:
        try:
            while True:
                time.sleep(1)
        except KeyboardInterrupt:
            pass


if __name__ == "__main__":
    main()
