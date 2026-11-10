"""Tests of the client against a fake server that speaks the folder protocol (no torch, no GPU needed).

    python -m unittest discover -s tests
"""

import io
import json
import os
import tempfile
import threading
import time
import unittest

import numpy as np
import tifffile

from lysoquant_client import LysoQuantError, check_server, client, make_rgb, segment


class FakeServer:
    """Writes the heartbeat and answers every job: labels of the size of the input, which is at the model pixel size"""

    def __init__(self, folder, fail=False, hang=False):
        self.folder, self.fail, self.hang = folder, fail, hang
        self.seen = []
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, daemon=True)

    def __enter__(self):
        self._beat()        # the server is up when the test starts
        self._thread.start()
        return self

    def __exit__(self, *exc):
        self._stop.set()
        self._thread.join()

    def _beat(self):
        path = os.path.join(self.folder, client.HEARTBEAT)
        with open(path + ".part", "w") as f:
            json.dump({"time": time.time(), "device": "cpu", "tile_size": 388, "folder": self.folder,
                       "job": None, "jobs": {}, "done": 0}, f)
        os.replace(path + ".part", path)

    def _run(self):
        while not self._stop.is_set():
            self._beat()
            for name in os.listdir(self.folder):
                job = os.path.join(self.folder, name)
                if name.startswith("lysoquant-") and name not in self.seen and os.path.exists(job + "/input.tif"):
                    self.seen.append(name)
                    self._answer(job)
            time.sleep(0.05)

    def _answer(self, job):
        if self.hang:
            return
        if self.fail:
            with open(job + "/error.txt", "w") as f:
                f.write("RuntimeError: out of memory\n")
            return
        with open(job + "/params.json") as f:
            params = json.load(f)
        image = tifffile.imread(job + "/input.tif")
        with open(job + "/progress.txt", "w") as f:
            f.write("1 1\n")
        h, w = image.shape[-3], image.shape[-2]
        self.params = params
        labels = np.full(((image.shape[0],) if image.ndim == 4 else ()) + (h, w), 2, np.uint16)
        el = params["element_size_um"][0]
        tifffile.imwrite(job + "/output.tif", labels, imagej=True, resolution=(1 / el, 1 / el),
                         metadata={"unit": "um"})
        with open(job + "/done.txt", "w") as f:
            f.write("0.1\n")


class ClientTest(unittest.TestCase):

    def setUp(self):
        self.folder = tempfile.mkdtemp()
        self.saved = (client.POLL_S, client.HEARTBEAT_WAIT_S)
        client.POLL_S, client.HEARTBEAT_WAIT_S = 0.1, 2.0

    def tearDown(self):
        client.POLL_S, client.HEARTBEAT_WAIT_S = self.saved

    def jobs_left(self):
        return [n for n in os.listdir(self.folder) if n.startswith("lysoquant-")]

    def test_segment(self):
        rgb = np.zeros((20, 30, 3), np.uint8)
        steps = []
        with FakeServer(self.folder):
            result = segment(rgb, 0.1, local=True, folder=self.folder, progress=lambda d, t: steps.append((d, t)))
        self.assertEqual(result.labels.shape, (80, 120))        # 0.1 um -> 0.025 um: 4 times as many pixels
        self.assertEqual(result.labels.dtype, np.uint16)
        self.assertAlmostEqual(result.pixel_size_um[0], 0.025, places=3)
        self.assertEqual(self.jobs_left(), [])

    def test_stack_and_pixel_size_pair(self):
        with FakeServer(self.folder):
            result = segment(np.zeros((3, 8, 8, 3), np.uint8), (0.2, 0.2), local=True, folder=self.folder)
        self.assertEqual(result.labels.shape, (3, 64, 64))

    def test_server_gets_the_image_at_the_model_pixel_size(self):
        with FakeServer(self.folder) as server:
            segment(np.zeros((10, 20, 3), np.uint8), (0.05, 0.1), local=True, folder=self.folder)
        self.assertEqual(server.params, {"element_size_um": [0.025, 0.025]})

    def test_finds_the_process_folder(self):
        with FakeServer(self.folder):
            info = check_server(local=True, folder=self.folder)
        self.assertEqual(info, {"device": "cpu", "tile_size": 388, "folder": self.folder})

    def test_server_error(self):
        with FakeServer(self.folder, fail=True):
            with self.assertRaisesRegex(LysoQuantError, "out of memory"):
                segment(np.zeros((8, 8, 3), np.uint8), 0.1, local=True, folder=self.folder)
        self.assertEqual(self.jobs_left(), [])

    def test_timeout_cancels_the_job(self):
        with FakeServer(self.folder, hang=True) as server:
            with self.assertRaises(TimeoutError):
                segment(np.zeros((8, 8, 3), np.uint8), 0.1, local=True, folder=self.folder, timeout=0.5)
        self.assertEqual(len(server.seen), 1)       # the job had been sent
        self.assertEqual(self.jobs_left(), [])      # and was deleted

    def test_cancel_marker_when_the_folder_cannot_be_deleted(self):
        class Stubborn(client.LocalConnection):
            def delete_folder(self, path):
                raise OSError("busy")
        client._connect_saved = client._connect
        client._connect = lambda *a: Stubborn()
        try:
            with FakeServer(self.folder, hang=True):
                with self.assertRaises(TimeoutError):
                    segment(np.zeros((8, 8, 3), np.uint8), 0.1, local=True, folder=self.folder, timeout=0.3)
        finally:
            client._connect = client._connect_saved
        (job,) = self.jobs_left()
        self.assertTrue(os.path.exists(os.path.join(self.folder, job, "cancel.txt")))   # tells the server to stop

    def test_no_server(self):
        with self.assertRaisesRegex(LysoQuantError, "not running"):
            segment(np.zeros((8, 8, 3), np.uint8), 0.1, local=True, folder=self.folder)

    def test_input_is_checked(self):
        with self.assertRaises(ValueError):
            segment(np.zeros((8, 8), np.uint8), 0.1, local=True, folder=self.folder)
        with self.assertRaises(ValueError):
            segment(np.zeros((8, 8, 3), np.float32), 0.1, local=True, folder=self.folder)


class RescaleTest(unittest.TestCase):

    def test_upsampling_is_bilinear_with_mirrored_edge(self):
        rgb = np.zeros((1, 2, 3), np.uint8)
        rgb[0, :, 0] = [0, 100]
        out = client.rescale_rgb(rgb, (0.05, 0.05))      # twice the pixels
        self.assertEqual(out.shape, (2, 4, 3))
        self.assertEqual(out[0, :, 0].tolist(), [0, 50, 100, 50])   # beyond the last pixel: mirrored
        self.assertEqual(out[1].tolist(), out[0].tolist())

    def test_stack_and_no_change_at_model_size(self):
        stack = np.random.default_rng(0).integers(0, 256, (2, 5, 6, 3), dtype=np.uint8)
        self.assertIs(client.rescale_rgb(stack, (0.025, 0.025)), stack)
        out = client.rescale_rgb(stack, (0.05, 0.075))
        self.assertEqual(out.shape, (2, 10, 18, 3))
        self.assertEqual(out[1].tolist(), client.rescale_rgb(stack[1], (0.05, 0.075)).tolist())


class MakeRgbTest(unittest.TestCase):

    def test_scaling(self):
        protein = np.array([[0, 100], [200, 400]], np.uint16)
        lysosome = np.array([[10, 10], [10, 10]], np.uint16)       # flat: black
        rgb = make_rgb(protein, lysosome)
        self.assertEqual(rgb.dtype, np.uint8)
        self.assertEqual(rgb.shape, (2, 2, 3))
        self.assertEqual(rgb[..., 0].tolist(), [[0, 64], [128, 255]])
        self.assertEqual(rgb[..., 1].max(), 0)
        self.assertEqual(rgb[..., 2].max(), 0)

    def test_display_range(self):
        protein = np.array([[0, 50, 100, 400]], np.uint16)
        rgb = make_rgb(protein, protein, protein_range=(0, 100), lysosome_range=(50, 100))
        self.assertEqual(rgb[0, :, 0].tolist(), [0, 128, 255, 255])
        self.assertEqual(rgb[0, :, 1].tolist(), [0, 0, 255, 255])

    def test_shapes_must_match(self):
        with self.assertRaises(ValueError):
            make_rgb(np.zeros((4, 4)), np.zeros((4, 5)))


if __name__ == "__main__":
    unittest.main()
