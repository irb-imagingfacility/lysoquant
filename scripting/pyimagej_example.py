"""LysoQuant from Python with pyimagej: segment and count the lysosomes of an image.

    pip install pyimagej
    python pyimagej_example.py image.tif

Fiji must have the LysoQuant plugin. Without a display (a cluster) the SSH key of the
server is needed (a password cannot be asked): pass key_file.
"""

import sys

import imagej

ij = imagej.init("sc.fiji:fiji", mode="headless")       # or the path of a Fiji.app

image = ij.IJ.openImage(sys.argv[1])

module = ij.command().run(
    "ch.irb.lysoquant.AnalyzeCommand",
    True,                                               # process: check the inputs
    {
        "image": image,
        "lysosomeChannel": 2,
        "proteinChannel": 3,
        "singleValues": True,
        # The server of Edit>Options>LysoQuant Settings... unless you name it:
        # "remote": False,                              # the server runs on this computer
        # "server": "gpu-server.example.org", "port": "22", "username": "lysoquant",
        # "keyFile": "/home/me/.ssh/lysoquant.key",
        # "processFolder": "/jobs",
    },
).get()

labels = module.getOutput("labels")                     # ImagePlus: 0 background, 1 empty, 2 loaded lysosomes
labels_array = ij.py.from_java(labels)                  # xarray / numpy
print("labels:", labels_array.shape, "pixel size (um):", labels.getCalibration().pixelWidth)

# The tables as pandas DataFrames, through CSV
import pandas as pd

module.getOutput("totals").saveAs("/tmp/lysoquant_totals.csv")
module.getOutput("singles").saveAs("/tmp/lysoquant_lysosomes.csv")
print(pd.read_csv("/tmp/lysoquant_totals.csv"))
print(pd.read_csv("/tmp/lysoquant_lysosomes.csv").groupby("Lysosome Type").size())

ij.dispose()
