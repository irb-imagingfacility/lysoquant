# LysoQuant

This is the binary of the LysoQuant plugin for Fiji.

This plugin requires a connection to a workstation (can be the local computer) running the
LysoQuant PyTorch segmentation server (see [server/README.md](../server/README.md)).
No other Fiji plugin is required.

## Installation

1. Set up the server on a workstation with a GPU and convert the model weights, as described in `server/README.md`.
1. Add LysoQuant plugin to Fiji, either through the update site or by placing the .jar file at {Fiji.app}/plugins
1. On ImageJ, Edit>Options>LysoQuant Settings... to configure the connection to the server, and check it with "Test connection" (the tile size is set on the server, see `server/README.md`). "Hide the segmentation image (LQ_)" there skips opening the segmentation window: the measurements are made as usual.

## Usage

1. Open a microscopy file, for example test.tif
1. Set the parameters to measure with Analyze> Set parameters...
1. Select a ROI corresponding to a cell to analyze. If no ROI is selected, all the image will be analyzed. ROIs can also be added to ROI Manager.
1. Run Analyze > LysoQuant. The user will be prompted with a selection for the channels. In this case, set 2 for the protein and 3 for the lysosomes
1. If the option for single values is unchecked, the image will be segmented and analyzed and summary values will be presented. If checked, also single values for each lysosome will be presented. In this case, single lysosomes will be measured with the specified range of channels and the measurements parameters defined above.

If you use this, please cite

**Deep learning approach for quantification of organelles and misfolded polypeptides delivery within degradative compartments**

Diego Morone, Alessandro Marazza, Timothy J. Bergmann, and Maurizio Molinari

_Molecular Biology of the Cell_ 2020

doi: [https://doi.org/10.1091/mbc.E20-04-0269](https://doi.org/10.1091/mbc.E20-04-0269)
