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

import org.scijava.ItemIO;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

/**
 * Segment the lysosomes of an image on the LysoQuant server. Without any dialog,
 * so it works from scripts and in headless Fiji. Example (Groovy):
 * <pre>
 * #@ ImagePlus image
 * #@output ImagePlus labels
 * labels = ij.command().run("ch.irb.lysoquant.SegmentCommand", true,
 *     "image", image, "lysosomeChannel", 2, "proteinChannel", 3).get().getOutput("labels")
 * </pre>
 *
 * @author Diego Morone
 */
@Plugin(type = org.scijava.command.Command.class, name = "LysoQuant Segment", headless = true)
public class SegmentCommand extends LysoQuantCommand {

    @Parameter(label = "Segmentation", type = ItemIO.OUTPUT,
               description = "16-bit labels: 0 background, 1 empty, 2 loaded lysosomes")
    ImagePlus labels;

    @Override
    public void run() {
        try {
            labels = analysis(new LysoQuantAnalysis()).segment(image, null, serverSettings());
        } catch (Exception e) {
            throw failure(e);
        }
    }
}
