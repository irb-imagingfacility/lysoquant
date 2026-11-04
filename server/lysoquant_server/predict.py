"""
Run LysoQuant U-Net segmentation on a job folder created by the Fiji client.

The job folder contains:
    input.tif     RGB image (stack), as sent by the LysoQuant plugin
    params.json   {"element_size_um": [y, x]}

While running the folder contains progress.txt ("<tiles done> <tiles in total>").
On success it will contain output.tif (16-bit label image at the model
resolution) and done.txt. On failure error.txt holds the message.
If the client cancels it removes the folder, after writing cancel.txt in it
(in case it cannot delete everything): the job is stopped at once, even in the
middle of a tile (between two layers of the network), and its folder deleted.

Usage:
    lysoquant-predict --weights lysoquant-7-16.pt [--gpu "GPU 0"] /path/to/cache/lysoquant-<uuid>/
"""

import argparse
import json
import math
import os
import shutil
import sys
import threading
import time
import traceback
import warnings
from pathlib import Path

import numpy as np
import tifffile
import torch

from .convert import load_model


PROGRESS_NAME = "progress.txt"      # "<tiles done> <tiles in total>", rewritten after each tile
CANCEL_NAME = "cancel.txt"          # written by the client when the user cancels
CANCELED = 2                        # run_job() result, besides 0 (done) and 1 (failed)
DEFAULT_TILE_SIZE = 900             # output tile size in px, see "Tile size and GPU memory" in README.md


class Canceled(Exception):
    """The client canceled the job: its folder is gone or holds cancel.txt"""


def check_canceled(folder):
    if not folder.is_dir() or (folder / CANCEL_NAME).exists():
        raise Canceled()


# The network is shared by the jobs of a GPU, so the cancel flag of the job
# is looked up per thread, in a hook called before each layer
_job = threading.local()


def _stop_if_canceled(module, inputs):
    event = getattr(_job, "canceled", None)
    if event is not None and event.is_set():
        raise Canceled()


def _watch_cancel(folder, canceled, finished):
    """Set canceled as soon as the client cancels; a tile takes too long to wait for its end."""
    while not finished.wait(0.2):
        try:
            check_canceled(folder)
        except Canceled:
            canceled.set()
            return
        except OSError:
            pass


def write_progress(folder, done, total):
    tmp = folder / f"{PROGRESS_NAME}.part"
    tmp.write_text(f"{done} {total}\n")
    tmp.replace(folder / PROGRESS_NAME)


def log(msg):
    name = threading.current_thread().name      # workers are named after their device
    print(msg if name == "MainThread" else f"[{name}] {msg}", flush=True)


# ----------------------------------------------------------------------
# Pre-processing. These replicate de.unifreiburg.unet.TrainingSample so
# that the network sees exactly the same data as with the caffe backend.

def to_frames(image):
    """Return a float32 array (frames, channels, y, x) from a TIFF array."""
    image = np.asarray(image)
    if image.ndim == 2:                                # single gray image
        image = image[None, None]
    elif image.ndim == 3 and image.shape[-1] in (3, 4):  # single RGB(A)
        image = np.moveaxis(image[..., :3], -1, 0)[None]
    elif image.ndim == 3:                              # gray stack
        image = image[:, None]
    elif image.ndim == 4 and image.shape[-1] in (3, 4):  # RGB(A) stack
        image = np.moveaxis(image[..., :3], -1, 1)
    else:
        raise ValueError(f"Unsupported image shape {image.shape}")
    return image.astype(np.float32)


def _interp_indices(n_in, n_out, scale):
    pos = np.arange(n_out) / scale
    lower = np.floor(pos).astype(np.int64)
    upper = lower + 1
    upper = np.where(upper < n_in, upper, 2 * (n_in - 1) - upper)
    return lower, upper, pos - lower


def rescale_xy(frames, scale_y, scale_x):
    """Bilinear rescaling with mirrored upper neighbours (TrainingSample.rescaleXY)."""
    if scale_y == 1 and scale_x == 1:
        return frames
    h, w = frames.shape[-2:]
    h_out, w_out = int(math.floor(h * scale_y + 0.5)), int(math.floor(w * scale_x + 0.5))
    yl, yu, dy = _interp_indices(h, h_out, scale_y)
    xl, xu, dx = _interp_indices(w, w_out, scale_x)
    dy = dy[:, None]
    f = frames.astype(np.float64)
    out = ((1 - dx) * (1 - dy) * f[..., yl[:, None], xl] +
           (1 - dx) * dy * f[..., yu[:, None], xl] +
           dx * (1 - dy) * f[..., yl[:, None], xu] +
           dx * dy * f[..., yu[:, None], xu])
    return out.astype(np.float32)


def normalize(frame, normalization_type):
    """Normalize one frame over all its channels (TrainingSample.normalizeValues)."""
    if normalization_type == 0:
        return frame, 1, 0
    if normalization_type == 1:      # [min, max] -> [0, 1]
        lo, hi = frame.min(), frame.max()
        offset = -float(lo)
        scale = 1.0 / float(hi - lo) if hi > lo else 1.0
    elif normalization_type == 2:    # zero mean, unit standard deviation
        offset = -float(frame.mean(dtype=np.float64))
        scale = math.sqrt(frame.size / np.sum((frame.astype(np.float64) + offset) ** 2))
    elif normalization_type == 3:    # max norm 1
        offset = 0.0
        scale = 1.0 / math.sqrt(np.max(np.sum(frame.astype(np.float64) ** 2, axis=0)))
    else:
        raise ValueError(f"Unknown normalization type {normalization_type}")
    offset, scale = np.float32(offset), np.float32(scale)
    return (frame + offset) * scale, scale, offset


# ----------------------------------------------------------------------
# Tiled prediction (caffe_unet tiled_predict with mirror padding)

def mirror_indices(start, stop, n):
    """Indices start..stop-1 mirrored (without edge repetition) into [0, n)."""
    idx = np.arange(start, stop)
    if n == 1:
        return np.zeros_like(idx)
    period = 2 * (n - 1)
    idx = np.mod(idx, period)
    return np.where(idx < n, idx, period - idx)


def tile_layout(size, requested, net):
    """Output tile size and number of tiles along one axis.

    When the whole axis fits into one tile we use the smallest valid tile
    that covers it: tile origins stay at 0, so the result is identical to
    using the bigger tile, just faster.
    """
    tile = net.valid_output_shape(requested)
    if size <= tile:
        tile = net.valid_output_shape(size)
        if tile < size:
            tile += net.downsample_factor
    return tile, math.ceil(size / tile)


@torch.inference_mode()
def predict_frame(net, frame, tile_size, device, progress="", on_tile=None):
    """Return class scores (classes, y, x) for one normalized frame (channels, y, x).

    on_tile(done, total) is called before the first tile and after each tile.
    """
    _, h, w = frame.shape
    tile_y, ny = tile_layout(h, tile_size, net)
    tile_x, nx = tile_layout(w, tile_size, net)
    pad_y = (net.input_shape(tile_y) - tile_y) // 2
    pad_x = (net.input_shape(tile_x) - tile_x) // 2

    # Mirror-pad the whole frame once so every tile can be sliced directly.
    # It stays in host memory: GPU memory only depends on the tile size.
    rows = mirror_indices(-pad_y, ny * tile_y + pad_y, h)
    cols = mirror_indices(-pad_x, nx * tile_x + pad_x, w)
    padded = torch.from_numpy(np.ascontiguousarray(frame[:, rows[:, None], cols]))

    scores = None
    if on_tile:
        on_tile(0, ny * nx)
    for iy in range(ny):
        for ix in range(nx):
            log(f"Segmenting {progress}tile {iy * nx + ix + 1}/{ny * nx}")
            y0, x0 = iy * tile_y, ix * tile_x
            tile = padded[None, :, y0:y0 + tile_y + 2 * pad_y, x0:x0 + tile_x + 2 * pad_x].to(device)
            out = net(tile)[0].float().cpu().numpy()
            if scores is None:
                scores = np.empty((out.shape[0], ny * tile_y, nx * tile_x), np.float32)
            scores[:, y0:y0 + tile_y, x0:x0 + tile_x] = out
            if on_tile:
                on_tile(iy * nx + ix + 1, ny * nx)
    return scores[:, :h, :w]


def select_device(gpu):
    gpu = (gpu or "none").strip()
    if gpu.lower() in ("none", "cpu", ""):
        return torch.device("cpu")
    if torch.cuda.is_available():
        index = gpu[4:] if gpu.upper().startswith("GPU ") else gpu if gpu.isdigit() else None
        if index is None:
            return torch.device("cuda")
        if int(index) >= torch.cuda.device_count():
            raise ValueError(f"'{gpu}' requested but only {torch.cuda.device_count()} CUDA device(s) found")
        return torch.device(f"cuda:{int(index)}")
    if torch.backends.mps.is_available():
        log(f"CUDA not available, using Apple MPS instead of '{gpu}'")
        return torch.device("mps")
    log(f"WARNING: no GPU available, running '{gpu}' request on CPU (slow)")
    return torch.device("cpu")


def select_devices(gpu):
    """Devices for a --gpu value: 'GPU 0', 'GPU 0,GPU 2' (or '0,2'), 'all' or 'none'."""
    gpu = (gpu or "none").strip()
    if gpu.lower() == "all":
        if torch.cuda.is_available():
            return [torch.device(f"cuda:{i}") for i in range(torch.cuda.device_count())]
        return [select_device("GPU 0")]
    devices = []
    for part in gpu.split(","):
        device = select_device(part)
        if device not in devices:           # e.g. several requests that all fall back to the CPU
            devices.append(device)
    return devices


def _load(weights, device):
    net, config = load_model(weights, device)
    net.apply(lambda module: module.register_forward_pre_hook(_stop_if_canceled))
    log(f"Loaded model '{config['name']}' on {device}")
    return net, config, device


def _skip_missing_nnpack():
    """Do not let PyTorch try NNPACK on CPUs without it (e.g. an emulated x86 VM without AVX2).

    The convolutions work without it (other CPU kernels are used), but each
    attempt prints "Could not initialize NNPACK! Reason: Unsupported hardware."
    """
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            available = torch._nnpack_available()
        if not available:
            torch.backends.nnpack.set_flags(False)
            log("NNPACK is not supported by this CPU, not used")
    except (AttributeError, RuntimeError):      # other PyTorch versions: keep its default
        pass


def _full_precision():
    _skip_missing_nnpack()
    torch.backends.cudnn.allow_tf32 = False      # keep full float32 precision
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.benchmark = True


def prepare_model(weights, gpu):
    """Load the model on the requested device; returns (net, config, device)."""
    _full_precision()
    return _load(weights, select_device(gpu))


def prepare_models(weights, gpu):
    """Load one copy of the model on each requested device; returns a list of prepare_model() results."""
    _full_precision()
    return [_load(weights, device) for device in select_devices(gpu)]


def segment(image, params, model, tile_size=DEFAULT_TILE_SIZE, progress=None):
    """Segment image with a prepare_model() result, in square tiles of tile_size px (output).

    progress(done, total) is called with the number of tiles done in all the frames.
    """
    net, config, device = model

    frames = to_frames(image)
    if frames.shape[1] != config["in_channels"]:
        raise ValueError(f"Model expects {config['in_channels']} channels, image has {frames.shape[1]}")

    el_data = [float(v) for v in params["element_size_um"]]
    el_model = config["element_size_um"]
    scale_y, scale_x = el_data[0] / el_model[0], el_data[1] / el_model[1]
    if scale_y != 1 or scale_x != 1:
        log(f"Rescaling Hyperstack (xy) from ({el_data[0]}, {el_data[1]}) to ({el_model[0]}, {el_model[1]})")

    labels = []
    for t, frame in enumerate(frames):
        frame = rescale_xy(frame, scale_y, scale_x)
        frame, scale, offset = normalize(frame, config["normalization_type"])
        log(f"t = {t + 1}: scale = {scale:g}, offset = {offset:g}")
        # all the frames have the same size, hence the same number of tiles
        on_tile = progress and (lambda done, total, t=t: progress(t * total + done, len(frames) * total))
        scores = predict_frame(net, frame, tile_size, device, f"frame {t + 1}/{len(frames)}, ", on_tile)
        labels.append(np.argmax(scores, axis=0).astype(np.uint16))  # first maximum wins, as in the U-Net plugin
    return np.stack(labels), el_model


def run_job(folder, model=None, weights=None, gpu=None, tile_size=DEFAULT_TILE_SIZE):
    """Segment a job folder with a loaded model, or load it from weights on gpu."""
    folder = Path(folder)
    start = time.time()
    canceled, finished = threading.Event(), threading.Event()
    _job.canceled = canceled
    threading.Thread(target=_watch_cancel, args=(folder, canceled, finished), daemon=True).start()
    try:
        check_canceled(folder)
        for name in ("done.txt", "error.txt", "output.tif", PROGRESS_NAME):
            (folder / name).unlink(missing_ok=True)
        params = json.loads((folder / "params.json").read_text())
        image = tifffile.imread(folder / "input.tif")
        def progress(done, total):      # called before the first tile and after each one
            check_canceled(folder)
            write_progress(folder, done, total)

        labels, el_model = segment(image, params, model or prepare_model(weights, gpu), tile_size, progress)

        tmp = folder / "output.tif.part"
        if labels.shape[0] == 1:
            labels = labels[0]
        tifffile.imwrite(tmp, labels, imagej=True,
                         resolution=(1.0 / el_model[1], 1.0 / el_model[0]),
                         metadata={"unit": "um", "axes": "TYX" if labels.ndim == 3 else "YX"})
        tmp.rename(folder / "output.tif")
        (folder / "done.txt").write_text(f"{time.time() - start:.1f}\n")
        log(f"Segmentation done in {time.time() - start:.1f} s")
        return 0
    except Canceled:
        shutil.rmtree(folder, ignore_errors=True)
        log(f"Job {folder.name} canceled by the client, stopped and deleted")
        return CANCELED
    except Exception as e:
        traceback.print_exc()
        (folder / "error.txt").write_text(f"{type(e).__name__}: {e}\n")
        return 1
    finally:
        finished.set()
        _job.canceled = None


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("folder", help="job folder containing input.tif and params.json")
    parser.add_argument("--weights", default=os.environ.get("LYSOQUANT_WEIGHTS"),
                        help="model weights (.pt) (default $LYSOQUANT_WEIGHTS)")
    parser.add_argument("--gpu", default=os.environ.get("LYSOQUANT_GPU", "GPU 0"),
                        help="'GPU <n>' or 'none' for CPU (default $LYSOQUANT_GPU or 'GPU 0')")
    parser.add_argument("--tile-size", type=int, default=int(os.environ.get("LYSOQUANT_TILE_SIZE", DEFAULT_TILE_SIZE)),
                        help=f"output tile size in px (default $LYSOQUANT_TILE_SIZE or {DEFAULT_TILE_SIZE})")
    args = parser.parse_args()
    if not args.weights:
        parser.error("--weights or $LYSOQUANT_WEIGHTS is required")
    sys.exit(run_job(args.folder, weights=args.weights, gpu=args.gpu, tile_size=args.tile_size))


if __name__ == "__main__":
    main()
