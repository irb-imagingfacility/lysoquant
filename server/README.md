# LysoQuant segmentation server (PyTorch)

PyTorch backend for the LysoQuant Fiji plugin. It replaces `caffe_unet` and
the U-Net Segmentation plugin; the Fiji side only needs the LysoQuant plugin.

The network is a one-to-one port of the original caffe_unet model, and the
pre-processing (bilinear rescaling to the model pixel size, min/max
normalization, mirror-padded tiling) replicates the U-Net plugin, so
segmentations match the Caffe backend.

## Deployments

`lysoquant-serve` keeps the model loaded, watches a process folder and segments
the jobs that the Fiji plugin puts there (see "Protocol"). Queued jobs run one
at a time per GPU by default, so with several GPUs (`--gpu all`) as many jobs run
at once as there are GPUs; `--tiles-per-gpu` (`$LYSOQUANT_TILES_PER_GPU`) raises
that, see "Tile size and GPU memory". Pick the one that fits:

| Where | How |
|---|---|
| Linux GPU server, users on other computers | [`../docker-ssh`](../docker-ssh/README.md): Docker image with an SFTP-only SSH server, access by public key |
| Windows PC (CUDA or CPU), one user or all users of the PC | [`../windows`](../windows/README.md) |
| Mac (Apple GPU or CPU) | [`../macos`](../macos/README.md) |
| Anything else | "Without Docker" below |

The server rewrites a heartbeat file, `.lysoquant-heartbeat.json`, every second
in the process folder, listing the running jobs. The Fiji plugin reads it and
reports an error at once if the server is down; `lysoquant-health` (used by
the Docker health check) reads it too: it fails if the heartbeat is older than
30 s or if a job runs longer than `LYSOQUANT_JOB_TIMEOUT`. Job folders left by
crashed clients are deleted after `LYSOQUANT_KEEP_HOURS` (default 24 h).

## Convert the Caffe model (once)

```bash
python -m venv ~/lysoquant-venv
~/lysoquant-venv/bin/pip install "/path/to/lysoquant-torch/server[convert]"
~/lysoquant-venv/bin/lysoquant-convert 2D-7-16.modeldef.h5 2D-7-16.caffemodel.h5 lysoquant-7-16.pt
```

Only `lysoquant-7-16.pt` is needed afterwards. The `.modeldef.h5` file is
no longer needed on the Fiji computer.

## Fiji settings (Edit > Options > LysoQuant Settings...)

| Field | Example |
|---|---|
| Server, port, username | `docker-ssh`: the host, its SSH port and `lysoquant` |
| Process folder on server | `docker-ssh`: `/jobs`; Windows and macOS: the process folder of the installer |
| SSH key | path to the private key, or empty to be asked for a password |

"Test connection" connects with the values in the dialog (without saving them),
checks that the server is running and shows the device and the tile size it uses.
Untick "Use remote host" if the server runs on the Fiji computer itself (Windows, macOS).

The tile size is a server setting, not a Fiji one (see below).

### Tile size and memory

Images are segmented in square tiles, whose size is set on the server with
`--tile-size` (`LYSOQUANT_TILE_SIZE` in the Docker `.env`; default 900 px).
Because of the unpadded convolutions, every tile needs 764 px of context on
each side (a 900 px tile reads 2428 px of input), so larger tiles are faster,
but they need more memory: choose the largest one that fits. Valid sizes are
4 + 128·k px (another value is rounded down to one of them).

Memory for one tile (float32, PyTorch): about 0.5 GB of weights, 0.4 GB of CUDA
runtime and the activations of the tile. cuDNN workspaces can add a few hundred
MB, so leave some margin. The image itself stays in host memory, so the GPU
memory does not depend on the image size. Without a GPU (`--gpu none`) the same
amount is needed in RAM.

| Tile (px) | Input read (px) | Memory (GB) |
|---|---|---|
| 132 | 1660 | 1.2 |
| 260 | 1788 | 1.3 |
| 388 | 1916 | 1.4 |
| 516 | 2044 | 1.4 |
| 644 | 2172 | 1.5 |
| 772 | 2300 | 1.6 |
| **900** (default) | 2428 | 1.7 |
| 1028 | 2556 | 1.8 |
| 1156 | 2684 | 1.9 |
| 1284 | 2812 | 2.1 |
| 1412 | 2940 | 2.2 |
| 1540 | 3068 | 2.4 |
| 1668 | 3196 | 2.6 |
| 1796 | 3324 | 2.8 |
| 1924 | 3452 | 3.0 |
| 2052 | 3580 | 3.3 |

By default each GPU runs one job at a time, so it only needs memory for one
tile. `--tiles-per-gpu N` (`LYSOQUANT_TILES_PER_GPU` in the Docker `.env`) lets
each GPU run N jobs at once, hence N tiles at a time (the tiles of one job are
still run one after the other): the GPU then needs memory for N tiles, but only
one copy of the weights and of the CUDA runtime (about 0.9 GB). Use it when the
GPU is big enough for several tiles and many users submit jobs together. If the
server runs out of memory, the job fails with a message suggesting a smaller
tile size.

## Without an installer

```bash
python -m venv ~/lysoquant-venv
~/lysoquant-venv/bin/pip install /path/to/lysoquant-torch/server
~/lysoquant-venv/bin/lysoquant-serve --weights lysoquant-7-16.pt ~/lysoquant/jobs
```

Install a CUDA build of PyTorch first if pip does not pick one
(see https://pytorch.org/get-started/locally/). `--gpu "GPU 1"` selects
another GPU, `--gpu none` the CPU. With several GPUs, `--gpu all` (or a list,
`--gpu "GPU 0,GPU 1"`) loads one copy of the model on each: queued jobs are
handed, oldest first, to whichever GPU is free, so as many jobs run at once as
there are GPUs. A single job still uses one GPU only.
A single job folder can also be segmented by hand with `lysoquant-predict --weights lysoquant-7-16.pt <job folder>`.

## Protocol

For each image the plugin checks the server heartbeat, creates the private
folder `<process folder>/lysoquant-<uuid>/` and writes `params.json`
(`element_size_um`) and `input.tif` (RGB: protein red,
lysosomes green), each as `.part` then renamed so that it appears complete.
The server picks the job up once `input.tif` exists, keeps `progress.txt`
(`<tiles done> <tiles in total>`) up to date for the plugin's progress bar and
writes `output.tif`
(16-bit labels: 0 background, 1 empty, 2 loaded) and `done.txt`, or
`error.txt` on failure. The plugin polls every second for these files,
then deletes the folder.
