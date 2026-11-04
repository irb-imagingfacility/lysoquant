/**************************************************************************
 *
 * Copyright (C) 2020 Diego Morone
 *
 *        Imaging Facility and Molinari Lab, 
 *	  Institute for Research in Biomedicine
 *	  Switzerland
 *	
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 **************************************************************************/


package ch.irb.lysoquant;

import ij.ImagePlus;
import ij.measure.ResultsTable;

import org.scijava.ItemIO;
import org.scijava.command.Command;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

/**
 * Count (and measure) the lysosomes of a segmentation made by LysoQuant
 * Segment. The slices and frames are the ones the segmentation was made with.
 * Needs no server.
 *
 * @author Diego Morone
 */
@Plugin(type = Command.class, name = "LysoQuant Count", headless = true)
public class CountCommand implements Command {

    @Parameter(label = "Image", description = "The image that was segmented")
    ImagePlus image;

    @Parameter(label = "Segmentation", description = "Result of LysoQuant Segment")
    ImagePlus labels;

    @Parameter(label = "Lysosome channel", min = "1")
    int lysosomeChannel = 2;

    @Parameter(label = "Protein channel", min = "1")
    int proteinChannel = 3;

    @Parameter(label = "First slice", min = "1")
    int firstSlice = 1;

    @Parameter(label = "Last slice (0: all)", min = "0")
    int lastSlice = 0;

    @Parameter(label = "First frame", min = "1")
    int firstFrame = 1;

    @Parameter(label = "Last frame (0: all)", min = "0")
    int lastFrame = 0;

    @Parameter(label = "First measurement channel", min = "1")
    int firstChannel = 1;

    @Parameter(label = "Last measurement channel (0: all)", min = "0")
    int lastChannel = 0;

    @Parameter(label = "Measure every lysosome")
    boolean singleValues = false;

    @Parameter(label = "Measure cargo area")
    boolean cargoArea = false;

    @Parameter(label = "Threshold for area", min = "0")
    int cargoThreshold = 0;

    @Parameter(label = "Filter min size", min = "0")
    double minSize = 0.53;

    @Parameter(label = "Totals", type = ItemIO.OUTPUT,
               description = "One row per image, slice and frame: number and fraction of lysosomes of each type")
    ResultsTable totals;

    @Parameter(label = "Lysosomes", type = ItemIO.OUTPUT,
               description = "One row per lysosome and measurement channel (empty without singleValues)")
    ResultsTable singles;

    @Override
    public void run() {
        LysoQuantAnalysis analysis = new LysoQuantAnalysis();
        analysis.lysoChannel = lysosomeChannel;
        analysis.proteinChannel = proteinChannel;
        analysis.firstC = firstChannel;
        analysis.lastC = lastChannel;
        analysis.firstZ = firstSlice;
        analysis.lastZ = lastSlice;
        analysis.firstT = firstFrame;
        analysis.lastT = lastFrame;
        analysis.singleValues = singleValues;
        analysis.cargoArea = cargoArea;
        analysis.cargoThreshold = cargoThreshold;
        analysis.minSize = minSize;
        analysis.check(image);

        // The segmentation of the selected slices and frames, in the order of a hyperstack
        int nz = analysis.lastZ - analysis.firstZ + 1;
        int nt = analysis.lastT - analysis.firstT + 1;
        if (labels.getStackSize() != nz * nt) {
            throw new IllegalArgumentException("The segmentation has " + labels.getStackSize() + " image(s), but "
                + nz + " slice(s) x " + nt + " frame(s) were selected");
        }
        labels.setDimensions(1, nz, nt);

        analysis.count(labels, image, null);
        totals = analysis.totals;
        singles = analysis.singles;
    }
}
