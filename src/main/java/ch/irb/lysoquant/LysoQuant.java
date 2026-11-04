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

import ij.ImageJ;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Roi;
import ij.gui.StackWindow;
import ij.measure.Measurements;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.plugin.filter.Analyzer;
import ij.plugin.frame.RoiManager;
import ij.text.TextWindow;
import ij.util.Tools;

import java.awt.Font;
import java.awt.Frame;
import java.awt.GraphicsEnvironment;
import java.awt.HeadlessException;

/**
 * LysoQuant deep learning segmentation of endolysosomes.
 * Takes 8-bit or 16-bit TIFF multichannel images (usually through Bioformats).
 * Deep learning model currently supports only 2D images
 * or stacks of uncorrelated distant slices.
 * <p>
 * The menu command: asks for the parameters, runs {@link LysoQuantAnalysis} and
 * shows the segmentation and the tables. Macros can fill in the dialog with
 * run("LysoQuant", "lysosome=2 protein=3 ..."). With the argument "segment"
 * (LysoQuant Segment...) it only segments, and shows the label image. Scripts can call
 * {@link LysoQuantAnalysis} or the SciJava commands directly.
 *
 * @author Diego Morone
 */
public class LysoQuant implements PlugIn, Measurements {
    int ch_lyso;
    int ch_protein;
    int firstC = 1;
    int lastC;
    int firstZ = 1;
    int lastZ;
    int nSlices;
    int firstT = 1;
    int lastT;
    int nFrames;
    int nChannels;
    boolean display_values;
    boolean display_cargo;
    int cargo_thres;
    boolean display_warning = Boolean.parseBoolean(Prefs.get("lysoquant.display_warning", "true"));
    boolean hide_segmentation = Boolean.parseBoolean(Prefs.get("lysoquant.hide_segmentation", "false"));

    double minSize = Double.parseDouble(Prefs.get("lysoquant.minsize", "0.53")); //in microns squared

    @Override
    public void run(String arg) {
        ImagePlus image = IJ.getImage();

        String problem = LysoQuantAnalysis.imageProblem(image);
        if (problem != null) {
            IJ.error(problem);
            return;
        }

        nChannels = image.getNChannels();
        lastC = nChannels;
        nSlices = image.getNSlices();
        lastZ = nSlices;
        nFrames = image.getNFrames();
        lastT = nFrames;

        boolean segmentOnly = "segment".equals(arg);
        if (!showDialog(segmentOnly)) return;

        // The cells of the ROI Manager, else the whole image (or its selection)
        Roi[] cells = null;
        RoiManager roiman = null;
        try {
            roiman = RoiManager.getInstance();
        } catch (HeadlessException e) {
            // no ROI Manager without a display (a macro in headless Fiji)
        }
        if (roiman != null && roiman.getCount() > 0) {
            roiman.runCommand(image, "Show None");
            cells = roiman.getRoisAsArray();
        }

        // Results are added to the tables that are open
        ResultsTable singles = openTable("Results");
        ResultsTable totals = openTable("LysoQuant");
        LysoQuantAnalysis analysis = new LysoQuantAnalysis(singles, totals);
        analysis.lysoChannel = ch_lyso;
        analysis.proteinChannel = ch_protein;
        analysis.firstC = firstC;
        analysis.lastC = lastC;
        analysis.firstZ = firstZ;
        analysis.lastZ = lastZ;
        analysis.firstT = firstT;
        analysis.lastT = lastT;
        analysis.singleValues = display_values;
        analysis.cargoArea = display_cargo;
        analysis.cargoThreshold = cargo_thres;
        analysis.minSize = minSize;
        analysis.addOverlay = true;

        ImagePlus segmented;
        try {
            segmented = analysis.segment(image, cells, ServerSettings.fromPrefs());
        } catch (InterruptedException e) {
            IJ.showStatus("LysoQuant: " + e.getMessage());
            return;
        } catch (Exception e) {
            IJ.error("LysoQuant", e.getMessage());
            return;
        } finally {
            IJ.showProgress(1.0);       // hide the progress bar
        }
        if (segmentOnly) {
            segmented.show();       // the result of the command, for the macro to select
            return;
        }
        if (!hide_segmentation) new StackWindow(segmented);

        analysis.count(segmented, image, cells);
        if (GraphicsEnvironment.isHeadless()) {
            // No windows: a macro reads the table with getResult(): the lysosomes, else the totals
            Analyzer.setResultsTable(display_values ? singles : totals);
            return;
        }
        if (display_values) singles.show("Results");
        totals.show("LysoQuant");
    }

    /** The table of the window with this title if it is open, else a new one */
    private static ResultsTable openTable(String title) {
        Frame frame = WindowManager.getFrame(title);
        if (frame instanceof TextWindow) {
            ResultsTable table = ((TextWindow) frame).getTextPanel().getResultsTable();
            if (table != null) return table;
        }
        return new ResultsTable();
    }

    /**
     * GUI for this plugin
     * 
     * @return true is everything goes fine, false if canceled
     */
    private boolean showDialog(boolean segmentOnly) {
        GenericDialog gd = new GenericDialog("LysoQuant");

	    gd.addNumericField("Lysosome Channel", Integer.parseInt(Prefs.get("lysoquant.display_lyso", "2")), 0);
	    gd.addNumericField("Protein Channel", Integer.parseInt(Prefs.get("lysoquant.display_protein", "3")), 0);

	    if (nSlices > 1) {
            gd.addStringField("Slices", "1-"+nSlices);
	    }	
	
	    if (nFrames > 1) {
            gd.addStringField("Frames", "1-"+nFrames);		      
        }	

	    if (!segmentOnly) {
	    boolean defaultTick_display = Boolean.parseBoolean(Prefs.get("lysoquant.display_values", "false"));
	    gd.addCheckbox("Show single values", defaultTick_display);

        int oldfirstC = (int)Prefs.get("lysoquant.display_firstC", 1);
        if (oldfirstC > nChannels) {
            oldfirstC = 1;
        }
        
        int oldlastC = (int)Prefs.get("lysoquant.display_lastC", nChannels);
        if (oldlastC > nChannels) {
            oldlastC = nChannels;
        }
    
        if (oldfirstC == oldlastC) {
            gd.addStringField("Measurement Channels", String.valueOf(oldfirstC));
        } else {
            gd.addStringField("Measurement Channels", String.valueOf(oldfirstC)+"-"+String.valueOf(oldlastC));
        }

        boolean defaultTick_area = Boolean.parseBoolean(Prefs.get("lysoquant.display_cargo", "false"));
        gd.addCheckbox("Measure cargo area", defaultTick_area);
        
        gd.addNumericField("Threshold for area", Integer.parseInt(Prefs.get("lysoquant.cargo_thres", "0")), 0); 
        }

        Font citationFont = new Font("Arial", Font.PLAIN, 10);
        gd.addMessage("Please cite Morone et al., MBoC 2020\ndoi:10.1091/mbc.e20-04-0269 ", citationFont);

        gd.showDialog();
        if (gd.wasCanceled())
            return false;

        // get entered values
        ch_lyso = (int)gd.getNextNumber();
        if (ch_lyso > nChannels) {
            IJ.error("Selected lysosome channel is higher than total number of channels");
            return false;
        }
        ch_protein = (int)gd.getNextNumber();
        if (ch_protein > nChannels) {
            IJ.error("Selected protein channel is higher than total number of channels");
            return false;
        }
        // The order of these ifs should match the order above
        if (nSlices > 1) {
            // Retrieve slices
            String[] range = Tools.split(gd.getNextString(), " -");
            double z1 = gd.parseDouble(range[0]);
            double z2 = range.length==2?gd.parseDouble(range[1]):Double.NaN;
            firstZ = Double.isNaN(z1)?1:(int)z1;
            lastZ = Double.isNaN(z2)?firstZ:(int)z2;
            if (firstZ<1) firstZ = 1;
            if (lastZ>nSlices) lastZ = nSlices;
            if (firstZ>lastZ) {firstZ=1; lastZ=nSlices;}
        }

        if (nFrames > 1) {
            // Retrieve frames
            String[] range = Tools.split(gd.getNextString(), " -");
            double t1 = gd.parseDouble(range[0]);
            double t2 = range.length==2?gd.parseDouble(range[1]):Double.NaN;
            firstT = Double.isNaN(t1)?1:(int)t1;
            lastT = Double.isNaN(t2)?firstT:(int)t2;
            if (firstT<1) firstT = 1;
            if (lastT>nFrames) lastT = nFrames;
            if (firstT>lastT) {firstT=1; lastT=nFrames;}
        }
	
        if (!segmentOnly) {
        display_values = gd.getNextBoolean();
        
        // Retrieve meas. channel range. Taken from Duplicator.java
        String[] range = Tools.split(gd.getNextString(), " -");
        double c1 = gd.parseDouble(range[0]);
        double c2 = range.length==2?gd.parseDouble(range[1]):Double.NaN;
        firstC = Double.isNaN(c1)?1:(int)c1;
        lastC = Double.isNaN(c2)?firstC:(int)c2;
        if (firstC<1) firstC = 1;
        if (lastC>nChannels) lastC = nChannels;
        if (firstC>lastC) {firstC=1; lastC=nChannels;}

        display_cargo = gd.getNextBoolean();

        cargo_thres = (int)gd.getNextNumber();
        }
       
        // Save for next usage
        Prefs.set("lysoquant.display_lyso", ch_lyso);
        Prefs.set("lysoquant.display_protein", ch_protein);
        if (!segmentOnly) {
            Prefs.set("lysoquant.display_values", Boolean.toString(display_values));
            Prefs.set("lysoquant.display_firstC", firstC);
            Prefs.set("lysoquant.display_lastC", lastC);
            Prefs.set("lysoquant.display_cargo", Boolean.toString(display_cargo));
            Prefs.set("lysoquant.cargo_thres", cargo_thres);
        }

        if (nSlices > 1 && display_warning && !IJ.isMacro()) {
            IJ.showMessage("LysoQuant is a 2D deep learning model", "Be careful! 3D images are not supported in this version of the deep learning model");
        }
        
        return true;
    }
    
    /**
     * Main method for debugging.
     *
     * For debugging, it is convenient to have a method that starts ImageJ, loads
     * an image and calls the plugin, e.g. after setting breakpoints.
     *
     * @param args unused
     */
    	public static void main(String[] args) {
		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = LysoQuant.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring("file:".length(), url.length() - clazz.getName().length() - ".class".length());
		System.setProperty("plugins.dir", pluginsDir);

		// start ImageJ
		new ImageJ();

		// open the Clown sample
		ImagePlus image = IJ.openImage("~/.lysoquant/mef.tif");
		image.show();

		// run the plugin
		IJ.runPlugIn(clazz.getName(), "");
	}
}
