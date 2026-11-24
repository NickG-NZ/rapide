import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import meshcat
import numpy as np
from view_recording import load_recording, matrix, populate, spring_pose


class RecordingTests(unittest.TestCase):
    def test_quaternion_order_and_translation(self):
        pose = [1, 2, 3, 0, 0, np.sqrt(.5), np.sqrt(.5)]
        np.testing.assert_allclose(matrix(pose) @ [1, 0, 0, 1], [1, 3, 3, 1], atol=1e-12)

    def test_rejects_irregular_samples(self):
        scene = {"version": 1, "fps": 50, "bodies": [], "frames": [{"time": .01, "poses": []}]}
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.json"
            path.write_text(json.dumps(scene))
            with self.assertRaises(ValueError):
                load_recording(path)

    @unittest.skipUnless(os.environ.get("RAPIDE_FOUR_BAR_EXECUTABLE"), "Set RAPIDE_FOUR_BAR_EXECUTABLE for C++ export integration")
    def test_export_to_meshcat_preserves_solver_poses(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "four_bar.json"
            subprocess.run([os.environ["RAPIDE_FOUR_BAR_EXECUTABLE"], "--record", str(path)], check=True, stdout=subprocess.DEVNULL)
            scene = load_recording(path)
        self.assertEqual(len(scene["frames"]), 501)
        self.assertAlmostEqual(scene["frames"][-1]["time"], 10)
        for frame in scene["frames"]:
            for joint in scene["joints"]:
                endpoints = []
                for side in ("a", "b"):
                    body = joint[side]
                    local = np.asarray(joint[f"point_{side}"]) - scene["bodies"][body]["com"]
                    endpoints.append((matrix(frame["poses"][body]) @ np.r_[local, 1])[:3])
                self.assertLess(np.linalg.norm(endpoints[0] - endpoints[1]), 2e-9)

        class Window:
            def __init__(self):
                self.commands = []

            def send(self, command):
                self.commands.append(command.lower())

        window = Window()
        populate(meshcat.Visualizer(window=window), scene)
        command = next(c for c in window.commands if c["type"] == "set_animation")
        self.assertEqual(len(command["animations"]), len(scene["bodies"]) + 3 * len(scene["springs"]))
        for body, animation in enumerate(command["animations"][:len(scene["bodies"])]):
            tracks = {track["name"]: track["keys"] for track in animation["clip"]["tracks"]}
            self.assertEqual(len(tracks[".position"]), 501)
            for index in (0, 250, 500):
                pose = scene["frames"][index]["poses"][body]
                np.testing.assert_allclose(tracks[".position"][index]["value"], pose[:3], atol=1e-12)
                np.testing.assert_allclose(matrix([*pose[:3], *tracks[".quaternion"][index]["value"]]), matrix(pose), atol=1e-12)

        spring = scene["springs"][0]
        coil_animation = next(a for a in command["animations"] if a["path"].endswith("/spring_0/coil"))
        scales = coil_animation["clip"]["tracks"][0]["keys"]
        lengths = []
        for index, frame in enumerate(scene["frames"]):
            transform, length = spring_pose(scene, frame, spring)
            lengths.append(length)
            np.testing.assert_allclose(scales[index]["value"], [1, length, 1], atol=1e-12)
            for side, local_y in (("a", 0), ("b", length)):
                body = spring[side]
                local = np.asarray(spring[f"point_{side}"]) - scene["bodies"][body]["com"]
                expected = matrix(frame["poses"][body]) @ np.r_[local, 1]
                np.testing.assert_allclose(transform @ [0, local_y, 0, 1], expected, atol=1e-12)
        self.assertGreater(max(lengths) - min(lengths), .05)


if __name__ == "__main__":
    unittest.main()
