"""LysoQuant cell by cell: Cellpose finds the cells, LysoQuant analyzes each one.

    pip install cellpose pyimagej tifffile pandas
    python cellpose_lysoquant.py [image.tif]

1. One channel, or several channels merged into an RGB, is segmented with a Cellpose model: one ROI per cell.
2. Every cell is cropped from the multichannel image, the rest of the crop is erased, and the crop is
   analyzed by LysoQuant (segmentation on the LysoQuant server, then counting), as in analyze.groovy.
3. The results of all the cells are collected in two tables (CSV) with a Cell column, next to the cells.

Fiji must have the LysoQuant plugin (FIJI_DIR), pyimagej needs Java 21 like Fiji (JAVA_HOME: the java folder
of Fiji works), and the LysoQuant server must be running (Edit>Options>LysoQuant
Settings... in Fiji, or LYSOQUANT_SERVER below). Without a display the SSH key of the server is needed
(a password cannot be asked).
"""

import sys
import tempfile
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy.ndimage import gaussian_filter

# ---------------------------------------------------------------------------------------------------
# Settings

IMAGE = Path(__file__).resolve().parent.parent / "test" / "TEST.tif"    # multichannel TIFF (c, y, x), or the first argument
OUTPUT_DIR = None               # default: <image name>_lysoquant next to the image

LYSOSOME_CHANNEL = 2            # channels count from 1, as in Fiji
PROTEIN_CHANNEL = 3

CELL_CHANNELS = [1, 2, 3, 4]    # the channels that show the cells: one channel (gray), or several merged into an RGB
                                # (each in the color of its LUT in the TIFF, as in Fiji's composite mode)

CELLPOSE_MODEL = "cpsam"        # name of a Cellpose model (cpsam, cyto3 with Cellpose 3, ...) or path to a custom one
DIAMETER = 300                  # of the cells in pixels (TEST.tif: one large cell); None: the model decides
SMOOTH_SIGMA_PX = 8             # blur the image for Cellpose: a cell full of lysosomes is not split into them. 0: no blur
FLOW_THRESHOLD = 0.4
CELLPROB_THRESHOLD = 0.0
USE_GPU = True                  # CUDA or Apple MPS if there is one
MIN_CELL_AREA_PX = 500          # smaller cells are not analyzed
EXCLUDE_EDGE_CELLS = True       # cells touching the border of the image are cut: not analyzed
CELL_PADDING_PX = 10            # margin around a cell in the crop sent to the server

FIJI_DIR = "/Applications/Fiji"             # the Fiji with the LysoQuant plugin (the folder with plugins/ and jars/)

# The server of Edit>Options>LysoQuant Settings... unless you name it here (the keys of AnalyzeCommand,
# see README.md): {"remote": False} for a server on this computer,
# {"server": "gpu-server.example.org", "port": "22", "username": "lysoquant", "keyFile": "/home/me/.ssh/lysoquant.key"}
LYSOQUANT_SERVER = {}

LYSOQUANT = {                   # the same as the options of Analyze>LysoQuant, see README.md
    "singleValues": True,       # measure every lysosome (the lysosomes table)
    "cargoArea": False,
    "cargoThreshold": 0,
    "minSize": 0.53,            # smallest lysosome, in squared microns
}
SAVE_LYSOSOME_LABELS = True     # the segmentation of every cell, as TIFF (0 background, 1 empty, 2 loaded lysosomes)

# ---------------------------------------------------------------------------------------------------

DEFAULT_COLORS = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1), (1, 1, 0)]    # when the TIFF has no LUTs


def read_channels(path):
    """The channels (c, y, x) of the TIFF, the color of each one (from its LUT), and its pixel size in microns"""
    with tifffile.TiffFile(path) as tif:
        series = tif.series[0]
        axes = "".join(a for a, n in zip(series.axes, series.shape) if n > 1)
        if axes != "CYX":
            raise ValueError(f"{path.name}: a multichannel image (c, y, x) is needed, found axes {series.axes} {series.shape}")
        channels = series.asarray().reshape([n for n in series.shape if n > 1])
        luts = (tif.imagej_metadata or {}).get("LUTs")
        if luts is not None and len(luts) == len(channels):
            colors = [tuple(np.asarray(lut)[:, -1] / 255.0) for lut in luts]
        else:
            colors = [DEFAULT_COLORS[i % len(DEFAULT_COLORS)] for i in range(len(channels))]
    return channels, colors


def cell_image(channels, colors, wanted):
    """What Cellpose sees: one channel (y, x), or the channels merged as an RGB (y, x, 3), as floats 0-1"""
    if any(not 1 <= c <= len(channels) for c in wanted):
        raise ValueError(f"CELL_CHANNELS {wanted}: the image has {len(channels)} channels")
    scaled = []
    for c in wanted:
        low, high = np.percentile(channels[c - 1], (1, 99.8))
        channel = np.clip((channels[c - 1].astype(np.float32) - low) / max(high - low, 1e-6), 0, 1)
        scaled.append(gaussian_filter(channel, SMOOTH_SIGMA_PX) if SMOOTH_SIGMA_PX else channel)
    if len(scaled) == 1:
        return scaled[0]
    rgb = sum(s[..., None] * np.asarray(colors[c - 1], np.float32) for s, c in zip(scaled, wanted))
    return np.clip(rgb, 0, 1).astype(np.float32)


def find_cells(image):
    """Label image of the cells (0 background, 1..n cells) with a Cellpose model"""
    from cellpose import models

    major = int(version("cellpose").split(".")[0])
    if major >= 4:
        model = models.CellposeModel(gpu=USE_GPU, pretrained_model=CELLPOSE_MODEL)
        options = {}
    else:       # Cellpose 3: models by name are model_type, custom ones are files; the channels are chosen above
        by_name = {"model_type": CELLPOSE_MODEL} if not Path(CELLPOSE_MODEL).exists() else {"pretrained_model": CELLPOSE_MODEL}
        model = models.CellposeModel(gpu=USE_GPU, **by_name)
        options = {"channels": [0, 0]}      # gray, or the mean of the three colors

    masks, _, _ = model.eval(image, diameter=DIAMETER, flow_threshold=FLOW_THRESHOLD,
                             cellprob_threshold=CELLPROB_THRESHOLD, **options)
    return masks.astype(np.uint16)


def select_cells(masks):
    """Labels of the cells to analyze, and why the others are left out"""
    border = set(np.unique(np.concatenate([masks[0], masks[-1], masks[:, 0], masks[:, -1]])))
    areas = np.bincount(masks.ravel())
    keep = []
    for label in range(1, len(areas)):
        if areas[label] == 0:
            continue
        if areas[label] < MIN_CELL_AREA_PX or (EXCLUDE_EDGE_CELLS and label in border):
            continue
        keep.append(label)
    return keep


def mask_to_roi(mask):
    """ImageJ selection (exact, holes included) of a boolean mask"""
    from jpype import JArray, JByte
    from scyjava import jimport

    ByteProcessor = jimport("ij.process.ByteProcessor")
    ThresholdToSelection = jimport("ij.plugin.filter.ThresholdToSelection")
    height, width = mask.shape
    pixels = (mask.astype(np.uint8) * 255).view(np.int8).ravel()
    ip = ByteProcessor(width, height, JArray(JByte)(pixels), None)
    ip.setThreshold(255, 255, ByteProcessor.NO_LUT_UPDATE)
    return ThresholdToSelection().convert(ip)


def read_table(table, path):
    """ResultsTable of ImageJ as a DataFrame (through CSV)"""
    table.saveAs(str(path))
    try:
        return pd.read_csv(path).drop(columns=[" "], errors="ignore")
    except pd.errors.EmptyDataError:        # no lysosome
        return pd.DataFrame()


def analyze_cells(path, masks, cells, output):
    """LysoQuant of every cell: the crop of the cell, with the rest erased, is segmented and counted"""
    import imagej
    from scyjava import jimport

    ij = imagej.init(FIJI_DIR, mode="headless")
    try:
        # The plugin erases what is outside of the cell with the background color of Fiji, white by default
        jimport("ij.gui.Toolbar").setBackgroundColor(jimport("java.awt.Color").BLACK)
        Duplicator = jimport("ij.plugin.Duplicator")

        image = ij.IJ.openImage(str(path))
        calibration = image.getCalibration()
        pixel_area = calibration.pixelWidth * calibration.pixelHeight
        height, width = masks.shape
        inputs_base = {
            "lysosomeChannel": LYSOSOME_CHANNEL,
            "proteinChannel": PROTEIN_CHANNEL,
            **LYSOQUANT,
            **LYSOQUANT_SERVER,
        }

        totals, lysosomes, failed = [], [], []
        scratch = Path(tempfile.mkdtemp(prefix="lysoquant_"))
        for n, label in enumerate(cells, 1):
            mask = masks == label
            ys, xs = np.nonzero(mask)
            x0, x1 = max(xs.min() - CELL_PADDING_PX, 0), min(xs.max() + 1 + CELL_PADDING_PX, width)
            y0, y1 = max(ys.min() - CELL_PADDING_PX, 0), min(ys.max() + 1 + CELL_PADDING_PX, height)

            image.setRoi(int(x0), int(y0), int(x1 - x0), int(y1 - y0))
            crop = Duplicator().run(image)
            image.deleteRoi()
            crop.setTitle(f"{path.stem}_cell{label:04d}")
            crop.setRoi(mask_to_roi(mask[y0:y1, x0:x1]))

            print(f"cell {label} ({n}/{len(cells)})", flush=True)
            try:
                module = ij.command().run("ch.irb.lysoquant.AnalyzeCommand", True, {"image": crop, **inputs_base}).get()
            except Exception as error:      # the server is not reachable: it would fail for every cell
                if not totals and not failed:
                    raise
                print(f"  failed: {error}", file=sys.stderr)
                failed.append(label)
                continue

            cell = {"Cell": label,
                    "Cell X (px)": round(float(xs.mean()), 1), "Cell Y (px)": round(float(ys.mean()), 1),
                    "Cell Area (um2)": round(len(xs) * pixel_area, 3)}
            for tables, name in ((totals, "totals"), (lysosomes, "singles")):
                frame = read_table(module.getOutput(name), scratch / f"{name}.csv")
                for i, (key, value) in enumerate(cell.items()):
                    frame.insert(i, key, value)
                tables.append(frame)
            if SAVE_LYSOSOME_LABELS:
                ij.IJ.saveAsTiff(module.getOutput("labels"), str(output / f"{path.stem}_cell{label:04d}_lysosomes.tif"))
    finally:
        ij.dispose()

    return pd.concat(totals, ignore_index=True), pd.concat(lysosomes, ignore_index=True), failed


def save_cells(output, stem, masks, shown):
    """Labels of the cells, and the image Cellpose saw with their outlines, to check the segmentation"""
    from cellpose import utils

    tifffile.imwrite(output / f"{stem}_cells.tif", masks)
    rgb = np.repeat(shown[..., None], 3, axis=-1) if shown.ndim == 2 else shown
    overlay = (rgb * 255).astype(np.uint8)
    overlay[utils.masks_to_outlines(masks)] = (255, 255, 255)
    tifffile.imwrite(output / f"{stem}_cells_overlay.tif", overlay, photometric="rgb")


def main():
    path = Path(sys.argv[1] if len(sys.argv) > 1 else IMAGE).resolve()
    output = Path(OUTPUT_DIR) if OUTPUT_DIR else path.parent / f"{path.stem}_lysoquant"
    output.mkdir(parents=True, exist_ok=True)

    channels, colors = read_channels(path)
    for name, channel in (("LYSOSOME_CHANNEL", LYSOSOME_CHANNEL), ("PROTEIN_CHANNEL", PROTEIN_CHANNEL)):
        if not 1 <= channel <= len(channels):
            raise ValueError(f"{name} {channel}: the image has {len(channels)} channels")

    shown = cell_image(channels, colors, CELL_CHANNELS)
    masks = find_cells(shown)
    save_cells(output, path.stem, masks, shown)
    cells = select_cells(masks)
    print(f"{masks.max()} cells found by Cellpose, {len(cells)} analyzed")
    if not cells:
        return

    totals, lysosomes, failed = analyze_cells(path, masks, cells, output)
    totals.to_csv(output / f"{path.stem}_cells_totals.csv", index=False)
    lysosomes.to_csv(output / f"{path.stem}_cells_lysosomes.csv", index=False)
    print(totals.to_string(index=False))
    if failed:
        print(f"Cells without result (see above): {failed}", file=sys.stderr)
    print(f"Results in {output}")


if __name__ == "__main__":
    main()
