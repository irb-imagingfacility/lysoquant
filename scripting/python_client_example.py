"""LysoQuant segmentation from Python without Fiji: see ../client/README.md

    pip install ../client
    python python_client_example.py image.tif
"""

import sys

import tifffile

from lysoquant_client import check_server, make_rgb, segment

server = dict(host="gpu-server.example.org", username="lysoquant", key_file="/home/me/.ssh/lysoquant.key")
# server = dict(local=True)                                     # the server runs on this computer

print(check_server(**server))                                   # device, tile size, process folder

stack = tifffile.imread(sys.argv[1])                            # channels (c, y, x)
rgb = make_rgb(protein=stack[2], lysosome=stack[1])             # channels 3 and 2: protein red, lysosomes green
result = segment(rgb, pixel_size_um=0.111, progress=lambda done, total: print(f"tile {done}/{total}"), **server)

labels = result.labels                                          # uint16: 0 background, 1 empty, 2 loaded
print(labels.shape, "pixel size (um):", result.pixel_size_um)
print("empty:", (labels == 1).sum(), "px, loaded:", (labels == 2).sum(), "px")
