# LysoQuant

This is the source of the LysoQuant plugin for Fiji.

This plugin requires the [U-Net Segmentation](https://github.com/lmb-freiburg/Unet-Segmentation) plugin and a connection to a Linux workstation (can be the local computer) running a special variant of caffe (caffe_unet).

The weight file for the U-Net network can be found in the releases

## Installation

1. Download the latest release.
1. Set up a workstation with a compiled version of caffe_unet (can be a docker file, as specified in the U-Net README).
1. Add the U-Net segmentation plugin and LysoQuant plugin to Fiji, through the update site.
1. Place the weight file on the workstation and the modelfile on the computer where Fiji is installed
1. Open the testRGB.tif image and run the Detection task in the U-Net segmentation plus. Set the parameters accordingly
1. Drag and drop setup.ijm, modify the paramters according to username and file locations, then run

## Usage

1. Open a microscopy file, for example test.tif
1. Set the parameters to measure with Analyze> Set parameters...
1. Select a ROI corresponding to a cell to analyze. If no ROI is selected, all the image will be analyzed
1. Run Analyze > LysoQuant. The user will be prompted with a selection for the channels. In this case, set 2 for the protein and 3 for the lysosomes
1. If the option for single values is unchecked, the image will be segmented and analyzed and summary values will be presented. If checked, also single values for each lysosome will be presented.

If you use this, please cite

**Deep learning approach for quantification of organelles and misfolded polypeptides delivery within degradative compartments**

Diego Morone, Alessandro Marazza, Timothy J. Bergmann, and Maurizio Molinari

_Molecular Biology of the Cell_ 2020

doi: [https://doi.org/10.1091/mbc.E20-04-0269](https://doi.org/10.1091/mbc.E20-04-0269)