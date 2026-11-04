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
 * Segment the lysosomes of an image on the LysoQuant server and count them:
 * LysoQuant Segment, then LysoQuant Count. Without any dialog, so it works from
 * scripts and in headless Fiji. Example (Jython):
 * <pre>
 * #@ ImagePlus image
 * #@output ImagePlus labels
 * #@output ResultsTable totals
 * </pre>
 * or from pyimagej: <code>ij.command().run("ch.irb.lysoquant.AnalyzeCommand", True, {...}).get()</code>.
 *
 * @author Diego Morone
 */
@Plugin(type = Command.class, name = "LysoQuant Analyze", headless = true)
public class AnalyzeCommand extends LysoQuantCommand {

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

    @Parameter(label = "Segmentation", type = ItemIO.OUTPUT,
               description = "16-bit labels: 0 background, 1 empty, 2 loaded lysosomes")
    ImagePlus labels;

    @Parameter(label = "Totals", type = ItemIO.OUTPUT,
               description = "One row per image, slice and frame: number and fraction of lysosomes of each type")
    ResultsTable totals;

    @Parameter(label = "Lysosomes", type = ItemIO.OUTPUT,
               description = "One row per lysosome and measurement channel (empty without singleValues)")
    ResultsTable singles;

    @Override
    public void run() {
        LysoQuantAnalysis analysis = analysis(new LysoQuantAnalysis());
        analysis.firstC = firstChannel;
        analysis.lastC = lastChannel;
        analysis.singleValues = singleValues;
        analysis.cargoArea = cargoArea;
        analysis.cargoThreshold = cargoThreshold;
        analysis.minSize = minSize;
        try {
            labels = analysis.analyze(image, null, serverSettings());
        } catch (Exception e) {
            throw failure(e);
        }
        totals = analysis.totals;
        singles = analysis.singles;
    }
}
