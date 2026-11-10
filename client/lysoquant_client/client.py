"""
Segment an image on a LysoQuant server, as the Fiji plugin does.

A job is a folder lysoquant-<uuid> in the process folder of the server: params.json
and input.tif are put in it, the server writes progress.txt while it works, then
output.tif and done.txt (or error.txt). Cancelling writes cancel.txt: the server
stops at once and deletes the folder.
"""

import io
import json
import math
import os
import time
import uuid
from typing import Callable, NamedTuple, Optional, Sequence, Union

import numpy as np
import tifffile

from .connection import Connection, LocalConnection, SFTPConnection

HEARTBEAT = ".lysoquant-heartbeat.json"
DEFAULT_FOLDERS = ("/jobs", "LysoQuant/jobs")       # where the server is looked for (relative to the SSH login)
MODEL_PIXEL_SIZE_UM = 0.025                         # pixel size of the model: images are rescaled to it here, not on the server
POLL_S = 1.0
HEARTBEAT_WAIT_S = 5.0


class LysoQuantError(RuntimeError):
    """The server is not running or the segmentation failed; the message is for the user"""


class Segmentation(NamedTuple):
    labels: np.ndarray                  # uint16 (y, x), or (n, y, x) for a stack: 0 background, 1 empty, 2 loaded
    pixel_size_um: "tuple[float, float]"  # (y, x) of the labels: the image is rescaled to the pixel size of the model


def make_rgb(protein, lysosome, protein_range=None, lysosome_range=None):
    """RGB image for segment(): protein in red, lysosomes in green.

    Each channel is scaled to 8 bit between its display range (the minimum and
    maximum of the channel by default), as Fiji does when it merges the channels.
    Give the display range you use in Fiji (min, max) to get the same input.

    protein, lysosome: arrays (y, x), or stacks (n, y, x) of independent images
    """
    protein, lysosome = np.asarray(protein), np.asarray(lysosome)
    if protein.shape != lysosome.shape:
        raise ValueError(f"The channels differ in shape: {protein.shape} and {lysosome.shape}")
    rgb = np.zeros(protein.shape + (3,), np.uint8)
    for i, (channel, lo_hi) in enumerate(((protein, protein_range), (lysosome, lysosome_range))):
        lo, hi = lo_hi if lo_hi is not None else (channel.min(), channel.max())
        scaled = (channel.astype(np.float64) - lo) / max(float(hi - lo), 1e-12)
        rgb[..., i] = np.round(np.clip(scaled, 0, 1) * 255)
    return rgb


def _interp_indices(n_in, n_out, scale):
    pos = np.arange(n_out) / scale
    lower = np.floor(pos).astype(np.int64)
    upper = lower + 1
    upper = np.where(upper < n_in, upper, 2 * (n_in - 1) - upper)
    return np.minimum(lower, n_in - 1), np.maximum(upper, 0), pos - lower


def rescale_rgb(rgb, pixel_size_um, target_um=MODEL_PIXEL_SIZE_UM):
    """Bilinear rescaling of an RGB image (y, x, 3), or a stack (n, y, x, 3), to the pixel size target_um.

    Pixel i of the result is at i / scale in the image and the neighbour beyond the last
    pixel is mirrored, as the Fiji plugin and the server do. The result is uint8.
    """
    scale_y, scale_x = pixel_size_um[0] / target_um, pixel_size_um[1] / target_um
    if scale_y == 1 and scale_x == 1:
        return rgb
    h, w = rgb.shape[-3:-1]
    h_out, w_out = int(math.floor(h * scale_y + 0.5)), int(math.floor(w * scale_x + 0.5))
    yl, yu, dy = _interp_indices(h, h_out, scale_y)
    xl, xu, dx = _interp_indices(w, w_out, scale_x)
    dy = dy[:, None, None]
    dx = dx[None, :, None]
    out = np.empty(rgb.shape[:-3] + (h_out, w_out, 3), np.uint8)
    for i in np.ndindex(*rgb.shape[:-3]):       # one image at a time: bounded memory
        f = rgb[i].astype(np.float64)
        top = f[yl[:, None], xl] * (1 - dx) + f[yl[:, None], xu] * dx
        bottom = f[yu[:, None], xl] * (1 - dx) + f[yu[:, None], xu] * dx
        out[i] = np.round(top * (1 - dy) + bottom * dy)
    return out


def _connect(host, port, username, key_file, password, local):
    if local or not host:
        return LocalConnection()
    return SFTPConnection(host, port, username, key_file, password)


def _folder_with_server(connection, folder):
    """The process folder: the one given, or the first usual one with a server in it"""
    candidates = [folder] if folder else list(DEFAULT_FOLDERS)
    if not folder and isinstance(connection, LocalConnection):
        candidates.append(os.path.join(os.path.expanduser("~"), "LysoQuant", "jobs"))
    for candidate in candidates:
        if connection.exists(_join(candidate, HEARTBEAT)):
            return candidate
    where = folder or "none of the usual process folders"
    raise LysoQuantError(f"The LysoQuant server is not running (no heartbeat in {where})")


def _join(folder, name):
    return folder.rstrip("/") + "/" + name


def _heartbeat(connection, folder):
    return json.loads(connection.get(_join(folder, HEARTBEAT)).decode("utf-8"))


def _check_server(connection, folder):
    """Fail early if the server is not running: busy, or its heartbeat changes within a few seconds"""
    down = f"The LysoQuant server is not running (no heartbeat in {_join(folder, HEARTBEAT)})"
    if not connection.exists(_join(folder, HEARTBEAT)):
        raise LysoQuantError(down)
    first = _heartbeat(connection, folder)
    if first.get("job") is not None:
        return first
    waited = 0.0
    while waited < HEARTBEAT_WAIT_S:
        time.sleep(POLL_S / 2)
        waited += POLL_S / 2
        latest = _heartbeat(connection, folder)
        if latest.get("time") != first.get("time"):
            return latest
    raise LysoQuantError(f"{down}: no heartbeat for {HEARTBEAT_WAIT_S:g} s")


def check_server(host=None, port=22, username=None, key_file=None, folder=None, password=None, local=False):
    """Connect and make sure the server is running.

    Returns what it reports about itself: {"device": ..., "tile_size": ..., "folder": ...}
    The process folder is looked for in the usual places if not given.
    """
    with _connect(host, port, username, key_file, password, local) as connection:
        folder = _folder_with_server(connection, folder)
        beat = _check_server(connection, folder)
        return {"device": beat.get("device"), "tile_size": beat.get("tile_size"),
                "folder": beat.get("folder", folder)}


def _tiff_bytes(array, **kwargs):
    buffer = io.BytesIO()
    tifffile.imwrite(buffer, array, **kwargs)
    return buffer.getvalue()


def _pixel_size_um(tiff_bytes):
    """(y, x) pixel size of an ImageJ TIFF written by the server (resolution in pixels per micron)"""
    with tifffile.TiffFile(io.BytesIO(tiff_bytes)) as tif:
        tags = tif.pages[0].tags
        sizes = []
        for name in ("YResolution", "XResolution"):
            num, den = tags[name].value
            sizes.append(den / num)
        return tuple(sizes)


def segment(
    rgb: np.ndarray,
    pixel_size_um: Union[float, Sequence[float]],
    host: Optional[str] = None,
    port: int = 22,
    username: Optional[str] = None,
    key_file: Optional[str] = None,
    folder: Optional[str] = None,
    password: Optional[str] = None,
    local: bool = False,
    timeout: Optional[float] = None,
    progress: Optional[Callable[[int, int], None]] = None,
) -> Segmentation:
    """Segment the lysosomes of an RGB image (protein in red, lysosomes in green) on the server.

    rgb:           uint8 array (y, x, 3), or a stack (n, y, x, 3) of independent images, e.g. from make_rgb()
    pixel_size_um: size of a pixel in microns: one number, or (y, x); the image is
                   rescaled to the model pixel size (0.025 um) here, so that the server does not have to
    host, port, username, key_file: the SSH server; key_file is the private key (no passphrase
                   asked). password is the alternative to a key. Without host, or with local=True,
                   the server is on this computer.
    folder:        process folder of the server; looked for in the usual places if not given
    timeout:       seconds to wait for the result, for the queue and the segmentation together
    progress:      called with (tiles done, tiles in total) while the server works

    Ctrl-C (or the timeout) cancels the job: the server stops at once and deletes it.
    Returns the labels (0 background, 1 empty, 2 loaded lysosomes) and their pixel size.
    """
    rgb = np.asarray(rgb)
    if rgb.dtype != np.uint8 or rgb.ndim not in (3, 4) or rgb.shape[-1] != 3:
        raise ValueError("rgb must be a uint8 array (y, x, 3) or (n, y, x, 3): see make_rgb()")
    if np.ndim(pixel_size_um) == 0:
        pixel_size_um = (pixel_size_um, pixel_size_um)
    size_y, size_x = (float(v) for v in pixel_size_um)
    if not (size_y > 0 and size_x > 0):
        raise ValueError("pixel_size_um must be positive")
    rgb = rescale_rgb(rgb, (size_y, size_x))
    params = json.dumps({"element_size_um": [MODEL_PIXEL_SIZE_UM] * 2}, indent=2).encode("utf-8") + b"\n"
    input_tiff = _tiff_bytes(rgb, photometric="rgb")

    with _connect(host, port, username, key_file, password, local) as connection:
        folder = _folder_with_server(connection, folder)
        _check_server(connection, folder)

        job = _join(folder, "lysoquant-" + str(uuid.uuid4()))
        connection.create_private_folder(job)
        canceled = True
        try:
            connection.put(params, job + "/params.json")
            connection.put(input_tiff, job + "/input.tif")      # the server takes the job once it is there

            deadline = None if timeout is None else time.monotonic() + timeout
            while not connection.exists(job + "/done.txt"):
                if connection.exists(job + "/error.txt"):
                    message = connection.get(job + "/error.txt").decode("utf-8", "replace").strip()
                    canceled = False
                    raise LysoQuantError("Segmentation failed on the server:\n" + message)
                if deadline is not None and time.monotonic() > deadline:
                    raise TimeoutError(f"No segmentation after {timeout:g} s")
                if progress is not None:
                    tiles = _read_progress(connection, job + "/progress.txt")
                    if tiles:
                        progress(*tiles)
                time.sleep(POLL_S)

            output = connection.get(job + "/output.tif")
            canceled = False
            return Segmentation(tifffile.imread(io.BytesIO(output)), _pixel_size_um(output))
        finally:
            if canceled:
                # The server may write in the folder while it is deleted and keep it alive:
                # the marker makes it stop at once and delete the folder itself
                try:
                    connection.put(b"", job + "/cancel.txt")
                except OSError:
                    pass
            try:
                connection.delete_folder(job)
            except OSError:
                pass            # the server deletes it, or the stale job cleanup does


def _read_progress(connection, path):
    try:
        if not connection.exists(path):
            return None
        done, total = (int(v) for v in connection.get(path).split()[:2])
        return (done, total) if total > 0 else None
    except (OSError, ValueError):
        return None             # being replaced, try again at the next poll
