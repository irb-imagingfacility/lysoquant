# lysoquant-client

Segment images with a LysoQuant server from Python, **without Fiji**. It speaks the same protocol as the
Fiji plugin: it puts the image in a job folder in the process folder of the server (over SSH/SFTP, or
on this computer), waits for the result, and deletes the folder. Counting and measuring the lysosomes is
done by the Fiji plugin (see [../scripting](../scripting/README.md)); this package returns the label image.

```
pip install ./client          # numpy, tifffile, paramiko
```

```python
import tifffile
from lysoquant_client import check_server, make_rgb, segment

server = dict(host="gpu-server.example.org", username="lysoquant", key_file="~/.ssh/lysoquant.key")
# server = dict(local=True)                  # the server runs on this computer

print(check_server(**server))                # {'device': 'cuda:0', 'tile_size': 900, 'folder': '/jobs'}

stack = tifffile.imread("image.tif")         # channels (c, y, x)
rgb = make_rgb(protein=stack[2], lysosome=stack[1])
result = segment(rgb, pixel_size_um=0.111, **server)
result.labels                                # uint16 (y, x): 0 background, 1 empty, 2 loaded lysosomes
result.pixel_size_um                         # (y, x) of the labels
```

- `segment(rgb, pixel_size_um, ...)`: `rgb` is uint8 `(y, x, 3)`, or a stack `(n, y, x, 3)` of independent
  images, with the protein in red and the lysosomes in green. The client rescales the image to the pixel
  size of the model (0.025 um) before sending it, so **the labels are larger or smaller than the image**;
  `result.pixel_size_um` is theirs. Rescaling also multiplies the data sent (16 times for 0.1 um images).
- `make_rgb(protein, lysosome)` scales each channel to 8 bit between its minimum and maximum, as
  Fiji does with the display range. For the same input as Fiji give the display range you use there:
  `protein_range=(min, max)`, `lysosome_range=(min, max)`.
- `folder` is the process folder of the server; it is looked for in the usual places (`/jobs`,
  `LysoQuant/jobs`) if not given. `key_file` is a private key without passphrase (or `password`).
- `timeout` (seconds) and Ctrl-C cancel the job: the server stops at once and deletes it.
- `progress=lambda done, total: ...` is called while the server works.
- Errors: `LysoQuantError` (server not running, segmentation failed on the server), `ConnectionError`
  (SSH), `TimeoutError`.

## Command line

```
lysoquant-segment image.tif -o labels.tif --pixel-size 0.111 --host gpu-server.example.org --username lysoquant --key ~/.ssh/lysoquant.key
lysoquant-segment --test --host gpu-server.example.org --username lysoquant --key ~/.ssh/lysoquant.key
```

Without `--host` the server is on this computer. The image is a TIFF with the channels as `(c, y, x)`
(`--lysosome 2 --protein 3`), or RGB with `--rgb`.

## Tests

```
python -m unittest discover -s tests         # against a fake server, no GPU needed
```
