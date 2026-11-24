import os
import csv
import io
from pathlib import Path
import subprocess
import tempfile
import unittest
from collections import Counter

import meshcat
import numpy as np
from tire_geometry import smooth_clearance, tire_mesh, tire_support
from meshcat import transformations
from view_recording import load_recording, populate


class TireGeometryTests(unittest.TestCase):
    def test_transition_is_smooth(self):
        blend = .02
        np.testing.assert_array_equal(smooth_clearance(np.linspace(-.1, 0, 20), blend), 0)
        self.assertAlmostEqual(float(smooth_clearance(blend, blend)), blend)
        step = 1e-6
        for boundary in (0, blend):
            left = smooth_clearance(boundary - np.arange(4) * step, blend)
            right = smooth_clearance(boundary + np.arange(4) * step, blend)
            self.assertAlmostEqual((3 * left[0] - 4 * left[1] + left[2]) / (2 * step),
                                   (-3 * right[0] + 4 * right[1] - right[2]) / (2 * step), delta=1e-5)
            self.assertAlmostEqual(np.dot([2, -5, 4, -1], left) / step**2,
                                   np.dot([2, -5, 4, -1], right) / step**2, delta=.01)

    def test_mesh_is_closed_and_outward(self):
        vertices, faces = tire_mesh(.35, .23, .14, .47)
        rings = vertices.reshape(96, 16, 3)
        radial = np.linalg.norm(rings[:, :, (0, 2)], axis=2)
        np.testing.assert_allclose(((radial - .29) / .06)**2 + (rings[:, :, 1] / .07)**2, 1, atol=1e-12)
        directed = Counter((int(a), int(b)) for face in faces for a, b in zip(face, np.roll(face, -1)))
        for (a, b), count in directed.items():
            self.assertEqual(count, 1)
            self.assertEqual(directed[b, a], 1)
        a, b, c = vertices[faces[:, 0]], vertices[faces[:, 1]], vertices[faces[:, 2]]
        self.assertGreater(np.sum(a * np.cross(b, c)) / 6, 0)

    def test_camber_and_tilted_plane_share_contact_geometry(self):
        for camber in (0, .4, -.6, np.pi / 2):
            pose = transformations.rotation_matrix(camber, [1, 0, 0])
            lowest = tire_support(.35, .23, .14, pose, np.array([0., 0., 1.]))
            pose[2, 3] = -lowest[2] - .015
            normal = np.array([0., 0., 1.])
            vertices, _ = tire_mesh(.35, .23, .14, pose=pose)
            self.assertGreaterEqual(vertices[:, 2].min(), -1e-12)
            self.assertLess(np.min(vertices[:, 2]), 1e-12)
            rotation = transformations.rotation_matrix(.6, [1, 2, 3])
            rotation[:3, 3] = [.2, -.4, .7]
            tilted = rotation @ pose
            new_normal = rotation[:3, :3] @ normal
            offset = np.dot(new_normal, rotation[:3, 3])
            rotated, _ = tire_mesh(.35, .23, .14, pose=tilted, ground_normal=new_normal, ground_height=offset)
            np.testing.assert_allclose(rotated, vertices @ rotation[:3, :3].T + rotation[:3, 3], atol=1e-12)
            support = tire_support(.35, .23, .14, tilted, new_normal)
            self.assertAlmostEqual(np.dot(support, new_normal) - offset, -.015)

    @unittest.skipUnless(os.environ.get("RAPIDE_WHEEL_EXECUTABLE"), "Set RAPIDE_WHEEL_EXECUTABLE for rollout integration")
    def test_recorded_tire_has_one_visible_shape_at_every_frame(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "wheel.json"
            result = subprocess.run([os.environ["RAPIDE_WHEEL_EXECUTABLE"], "--camber", ".35", "--record", str(path)], check=True, capture_output=True, text=True)
            scene = load_recording(path)
        self.assertEqual(len(scene["frames"]), 301)
        self.assertEqual(len(scene["tires"]), 1)
        from view_recording import matrix
        samples = list(csv.DictReader(io.StringIO(result.stdout)))
        tire = scene["tires"][0]
        normal = np.asarray(tire["ground_normal"])
        for index, frame in enumerate(scene["frames"][1:], start=1):
            pose = matrix(frame["poses"][tire["body"]])
            lowest = tire_support(tire["radius"], tire["rim_radius"], tire["width"], pose, normal)
            gap = np.dot(lowest, normal) - tire["ground_height"]
            sample = samples[2 * index - 1]
            np.testing.assert_allclose(lowest - gap * normal,
                [float(sample[f"patch_{axis}"]) for axis in "xyz"], atol=1e-9)
            self.assertAlmostEqual(max(0, -gap), float(sample["compression"]), delta=1e-9)


        class Window:
            def __init__(self):
                self.commands = []

            def send(self, command):
                self.commands.append(command.lower())

        window = Window()
        populate(meshcat.Visualizer(window=window), scene)
        command = next(c for c in window.commands if c["type"] == "set_animation")
        samples = [a for a in command["animations"] if "/tire_0/sample_" in a["path"]]
        self.assertEqual(len(samples), 301)
        for index in range(301):
            visible = []
            for sample in samples:
                keys = sample["clip"]["tracks"][0]["keys"]
                value = [key["value"] for key in keys if key["time"] <= index][-1]
                if value:
                    visible.append(int(sample["path"].split("sample_")[-1]))
            self.assertEqual(visible, [index])


if __name__ == "__main__":
    unittest.main()
