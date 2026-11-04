# Scripting LysoQuant

LysoQuant can be run without its dialogs from ImageJ macros, Fiji scripts (Groovy, Jython, ...),
Python (pyimagej), headless Fiji, and from Python without Fiji. All of them use the same segmentation
server as the menu command: **Edit>Options>LysoQuant Settings...** (or the parameters below) says where it is.

| From | Use | Example |
|---|---|---|
| ImageJ macro | `run("LysoQuant", "...")`, `run("LysoQuant Segment...", "...")` | [macro_example.ijm](macro_example.ijm) |
| Fiji script (Groovy, Jython) | the commands `LysoQuant Analyze`, `LysoQuant Segment`, `LysoQuant Count` | [analyze.groovy](analyze.groovy) |
| Python with Fiji | the same commands through pyimagej | [pyimagej_example.py](pyimagej_example.py) |
| Python: Cellpose cells, then LysoQuant in each | pyimagej, `LysoQuant Analyze` per cell | [cellpose_lysoquant.py](cellpose_lysoquant.py) |
| Python without Fiji | the `lysoquant-client` package, [../client](../client/README.md) | [python_client_example.py](python_client_example.py) |

## Macros

`run("LysoQuant", "lysosome=2 protein=3 slices=1-5 frames=1-3 show measurement=1-4 measure threshold=100")`
fills in the dialog of **Analyze>LysoQuant**. The keys are the first words of its fields (a checkbox is
ticked if its key is there): `lysosome`, `protein`, `slices` and `frames` (only for stacks and time
series), `show` (single values), `measurement` (channels), `measure` (cargo area), `threshold`.
The Macro Recorder records the call. The results are the same as with the dialog: the `Results` and
`LysoQuant` tables, and the `LQ_` image unless it is hidden in the settings.

`run("LysoQuant Segment...", "lysosome=2 protein=3")` only segments: the new image is the labels
(16-bit: 0 background, 1 empty, 2 loaded lysosomes). The cells of the ROI Manager, or the selection of
the image, restrict the segmentation, as with the menu command.

The settings are in `Prefs`: `call("ij.Prefs.set", "lysoquant.server", "host")`, see the example.

## Commands (Groovy, Jython, pyimagej)

Three SciJava commands, run by their class name with `CommandService.run(name, true, inputs...)`.
They open no dialog and no window, so they work in headless Fiji. Outputs are read from the module.

| Command | Does | Outputs |
|---|---|---|
| `ch.irb.lysoquant.AnalyzeCommand` | segments, then counts | `labels` (ImagePlus), `totals`, `singles` (ResultsTables) |
| `ch.irb.lysoquant.SegmentCommand` | segments | `labels` |
| `ch.irb.lysoquant.CountCommand` | counts a segmentation, no server needed | `totals`, `singles` |

Inputs (the names to use in the inputs):

- `image`: the multichannel image (8 or 16 bit); its selection restricts the analysis.
- `lysosomeChannel`, `proteinChannel`: from 1 (defaults 2 and 3).
- `firstSlice`, `lastSlice`, `firstFrame`, `lastFrame`: from 1; 0 as last is up to the end.
- Analyze and Count: `firstChannel`, `lastChannel` (channels to measure), `singleValues` (measure every
  lysosome: the `singles` table), `cargoArea`, `cargoThreshold`, `minSize` (microns squared, default 0.53).
- Count: `labels`, made by Segment with the same slices and frames.
- Analyze and Segment: the server, if not the one of the settings: `remote` (false: this computer),
  `server`, `port`, `username`, `keyFile`, `processFolder`.

`totals` has one row per image, slice and frame (counts and fractions of empty and loaded lysosomes);
`singles` one row per lysosome and measurement channel. **The measurements in `singles` (Area, Mean, ...)
are the ones ticked in Analyze>Set Measurements...**, as with the menu command: in headless Fiji the
settings of the Fiji installation apply.

## Cell by cell with Cellpose

[cellpose_lysoquant.py](cellpose_lysoquant.py) takes a multichannel TIFF, finds the cells with a Cellpose model
(one channel, or several merged into an RGB: `CELL_CHANNELS`), and runs `AnalyzeCommand` on every cell: the cell is
cropped, the rest of the crop is erased (black), and the lysosome and protein channels (`LYSOSOME_CHANNEL`,
`PROTEIN_CHANNEL`, at the top of the script) go to the server. One server job per cell. The output folder has
the cell labels and an overlay to check them, `*_cells_totals.csv` (one row per cell) and `*_cells_lysosomes.csv`
(one row per lysosome and measurement channel), both with the `Cell` number, position and area, and optionally
the lysosome segmentation of every cell.

```
pip install cellpose pyimagej tifffile pandas
export JAVA_HOME=/Applications/Fiji/java/macos-arm64/zulu21...     # the Java of Fiji, which current Fiji plugins need
python cellpose_lysoquant.py image.tif
```

The `DIAMETER` and `SMOOTH_SIGMA_PX` defaults suit `test/TEST.tif` (one large cell full of lysosomes). For
images with many small cells, set the diameter of a cell in pixels (or `None`) and no blur: Cellpose splits a
textured cell into pieces that look like cells itself. Look at `*_cells_overlay.tif` before trusting the counts.

## Headless

- The commands need no display. A **password cannot be asked**: give the SSH key (`keyFile` or the
  settings; without passphrase), or use a server on the same computer (`remote` false).
- A macro in headless Fiji works too, but shows nothing: the `Results` table holds the lysosomes (with
  `show`, else the totals) for `getResult()`; there is no `LysoQuant` table. Scripts are easier.

## Cancelling

Esc in Fiji stops the segmentation. A script is canceled by interrupting its thread (stopping the
script, Ctrl-C in pyimagej). The server stops the job at once and deletes its folder.

## Notes

- Only the part of the stack selected (slices and frames) is sent to the server and segmented.
- The label image has the pixel size of the model (not the one of the image), with the calibration set.
- Commands, macros and the menu command give the same counts: they share one implementation
  (`LysoQuantAnalysis`, usable directly from Java).
