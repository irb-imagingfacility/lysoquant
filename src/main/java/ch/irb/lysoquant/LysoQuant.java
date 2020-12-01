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

import java.awt.*;
import ij.*;
import ij.text.*;
import ij.gui.GenericDialog;
import ij.gui.Roi;
import ij.Prefs;
import ij.plugin.ChannelSplitter;
import ij.plugin.PlugIn;
import ij.plugin.RoiScaler;
import ij.plugin.RGBStackConverter;
import ij.plugin.RGBStackMerge;
import de.unifreiburg.unet.*;
import ij.measure.Calibration;
import ij.measure.Measurements;
import ij.measure.ResultsTable;
import ij.plugin.filter.Analyzer;
import ij.plugin.filter.ParticleAnalyzer;
import ij.plugin.frame.RoiManager;
import ij.process.ImageProcessor;
import static java.lang.Math.floor;
import ij.util.Tools;
import ij.plugin.HyperStackConverter;

import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * LysoQuant deep learning segmentation of endolysosomes.
 * Takes 8-bit or 16-bit TIFF multichannel images (usually through Bioformats).
 * Deep learning model currently supports only 2D images
 * or stacks of uncorrelated distant slices.
 *
 * @author Diego Morone
 */
public class LysoQuant implements PlugIn, Measurements {
    int ch_lyso;
    int ch_protein;
    int firstC = 1;
    int lastC;
    int nChannels;
    int firstZ = 1;
    int lastZ;
    int nSlices;
    int firstT = 1;
    int lastT;
    int nFrames;
    boolean display_values;
    boolean display_warning = Boolean.parseBoolean(Prefs.get("lysoquant.display_warning", "true"));
    String positives = "loaded";
    int pos = 2;
    String negatives = "empty";
    int neg = 1;
    double minSize = Double.parseDouble(Prefs.get("lysoquant.minsize", "")); //in microns squared
    String modelname = Prefs.get("lysoquant.model", "");

    String modelpath = Prefs.get("unet.modelDefinitionFolder", "")+"/"+modelname;
    String tilesize = Prefs.get("lysoquant.tilesize", "");
    String weightspath = Prefs.get("lysoquant.weights", "");
    String gpuflag = Prefs.get("lysoquant.gpuflag", "");
    String useremotehost = Prefs.get("lysoquant.remotehost", "");
    String hostname = Prefs.get("lysoquant.server", "");
    String port = Prefs.get("lysoquant.port", "");
    String username = Prefs.get("lysoquant.username", "");
    String keypath = Prefs.get("unet.rsaKeyFilename", "");
    String cachefolder = Prefs.get("unet.processfolder", "");
    String averageflag = "none";
    String keeporiginal = "false";
    String outputscores = "false";
    String outputsoftmaxscores = "false";
        
    @Override
    public void run(String arg) {

        ImagePlus image = IJ.getImage();

        CompositeImage cimg = image.isComposite()?(CompositeImage)image:null;
        int type = image.getType();
        if (cimg==null && !(type==ImagePlus.GRAY8 || type==ImagePlus.GRAY16)) {
            IJ.error("8-bit or 16-bit grayscale stack required");
            return;
        }
        int size = image.getStackSize();
        if (size<2 && cimg==null) {
            IJ.error("A minimum of 2 channels are required");
            return;
        }

        nChannels = image.getNChannels();
        lastC = nChannels;

        nSlices = image.getNSlices();
        lastZ = nSlices;

        nFrames = image.getNFrames();
        lastT = nFrames;

        Roi roiA = image.getRoi();
        RoiManager roiman = RoiManager.getInstance();

        String title = image.getTitle();

        if (showDialog()) {
            // parameters
            String unetp = "";
            unetp += "modelFilename=";
            unetp += modelpath;
            unetp += ",Tile shape (px):=";
            unetp += tilesize+"x"+tilesize;
            unetp += ",weightsFilename=";
            unetp += weightspath;
            unetp += ",gpuId=";
            unetp += gpuflag;
            unetp += ",useRemoteHost=";
            unetp += useremotehost;
            unetp += ",hostname=";
            unetp += hostname;
            unetp += ",port=";
            unetp += port;
            unetp += ",username=";
            unetp += username;
            unetp += ",RSAKeyFile=";
            unetp += keypath;
            unetp += ",processFolder=";
            unetp += cachefolder;
            unetp += ",average=";
            unetp += averageflag;
            unetp += ",keepOriginal=";
            unetp += keeporiginal;
            unetp += ",outputScores=";
            unetp += outputscores;
            unetp += ",outputSoftmaxScores=";
            unetp += outputsoftmaxscores;

            ImagePlus rgb = make_rgb(ch_protein, ch_lyso, image);
            rgb.show();

            if (roiman==null || roiman.getCount()<1){

                if (roiA != null) {
                    rgb.setRoi(roiA);
                    IJ.run("Clear Outside");
                }
                
                try {
                    SegmentationJob.processHyperStack(unetp);
                } catch (InterruptedException ex) {
                    Logger.getLogger(LysoQuant.class.getName()).log(Level.SEVERE, null, ex);
                }

                // FIXME: hyperstack conversion
                ImagePlus segmented = IJ.getImage();
                if (nFrames > 1 || nSlices > 1)
                    segmented = HyperStackConverter.toHyperStack(segmented, 1, nSlices, nFrames);

                segmented.setTitle("LQ_"+title);

                int totalpos = count(segmented, image, firstC, lastC, pos, positives, minSize, display_values);
                int totalneg = count(segmented, image, firstC, lastC, neg, negatives, minSize, display_values);

                ResultsTable totals = null;
                Frame frame = WindowManager.getFrame("LysoQuant");
                if (frame!=null && (frame instanceof TextWindow)) {
                TextWindow tw = (TextWindow)frame;
                ResultsTable table = tw.getTextPanel().getResultsTable();
                    if (table!= null) {
                        totals = table;
                    } else {
                        totals = new ResultsTable();
                    }
                } else {
                    totals = new ResultsTable();
                }

                totals.incrementCounter();
                totals.addLabel(title);
                totals.addValue("Loaded", totalpos);
                totals.addValue("Empty", totalneg);
                totals.addValue("Ratio", (double)totalpos/(double)(totalpos+totalneg));
                totals.addValue("Lysosome Ch", ch_lyso);
                totals.addValue("Protein Ch", ch_protein);
                totals.show("LysoQuant");

            } else {
                roiman.runCommand("Combine");
                IJ.run("Clear Outside");
                
                try {
                    SegmentationJob.processHyperStack(unetp);
                } catch (InterruptedException ex) {
                    Logger.getLogger(LysoQuant.class.getName()).log(Level.SEVERE, null, ex);
                }

                // FIXME: hyperstack conversion
                ImagePlus segmented = IJ.getImage();
                if (nFrames > 1 || nSlices > 1)
                    segmented = HyperStackConverter.toHyperStack(segmented, 1, nSlices, nFrames);

                segmented.setTitle("LQ_"+title);

                int width = image.getWidth();
                int swidth = segmented.getWidth();
                double scale = (double) swidth/ (double) width;

                Roi[] rois = roiman.getRoisAsArray();
                for (Roi roi : rois){
                    int pan_x = roi.getBounds().x;
                    int pan_y = roi.getBounds().y;
                        
                    int new_x = (int) floor(pan_x * scale);
                    int new_y = (int) floor(pan_y * scale);
                        
                    Roi scaled = RoiScaler.scale(roi, scale, scale, false);
                    scaled.setLocation(new_x, new_y);
                    String roiname = roi.getName();

                    int totalpos = count(segmented, image, scaled, firstC, lastC, pos, positives, minSize, display_values);
                    int totalneg = count(segmented, image, scaled, firstC, lastC, neg, negatives, minSize, display_values);
        
                    ResultsTable totals = null;
                    Frame frame = WindowManager.getFrame("LysoQuant");
                    if (frame!=null && (frame instanceof TextWindow)) {
                        TextWindow tw = (TextWindow)frame;
                        ResultsTable table = tw.getTextPanel().getResultsTable();
                        if (table!= null) {
                            totals = table;
                        } else {
                            totals = new ResultsTable();
                        }
                    } else {
                        totals = new ResultsTable();
                    }
                
                    totals.incrementCounter();
                    totals.addLabel(title+"-"+roiname);
                    totals.addValue("Loaded", totalpos);
                    totals.addValue("Empty", totalneg);
                    totals.addValue("Ratio", (double)totalpos/(double)(totalpos+totalneg));
                    totals.addValue("Lysosome Ch", ch_lyso);
                    totals.addValue("Protein Ch", ch_protein);
                    totals.show("LysoQuant");
                }
            }
        }
    }
    
    /**
     * Take the segmented image and count the number or objects in selected objClass
     * If display_values is true, also measure on each object depending on what is selected
     * in Analyze>Set measurements...
     * 
     *  @param segmented image to extract the ROIs
     *  @param raw image to measure
     *  @param firstC measurement channel
     *  @param lastC measurement channel
     *  @param objClass the values of output in segmented images. See class variables
     *  @param objName their corresponding names. See class variables
     *  @param minSize cutoff for recognizing particles
     *  @param display_values if true measure the value on each lysosome
     *  @return int of total count for objClass
     *  @return table with measurements and info about inputs and image position in hyperstack
     */
    private int count(ImagePlus segmented, ImagePlus raw, int firstC, int lastC, 
                        int objClass, String objName, double minSize, boolean display_values) {
        int width = raw.getWidth();
        int swidth = segmented.getWidth();
        double invscale = (double) width/ (double) swidth;

        int total = 0;

        // Constructors
        RoiManager countman = new RoiManager(true); // Hidden roimanager for this task
        ResultsTable singles = null; // Table for single results
        Frame frame = WindowManager.getFrame("Results");
        if (frame!=null && (frame instanceof TextWindow)) {
            TextWindow tw = (TextWindow)frame;
            ResultsTable table = tw.getTextPanel().getResultsTable();
            if (table!= null) {
                singles = table;
            }
        } else {
                singles = new ResultsTable();
        }
        Analyzer measure = new Analyzer(raw, singles);

        // Get segmented image and apply binary threshold for positives or negatives
        ImageProcessor ip = segmented.getProcessor();
        Calibration cal = segmented.getCalibration();
        ip.setThreshold((double) objClass, (double) objClass, ImageProcessor.NO_LUT_UPDATE);

        // Get the Rois from the segmented image
        double unitSquared = cal.pixelWidth*cal.pixelHeight;
        minSize = minSize / unitSquared; // minsize must be a double in pixel units
        int options = 0;
        options += ParticleAnalyzer.SHOW_NONE;
        options += ParticleAnalyzer.ADD_TO_MANAGER;
        options += ParticleAnalyzer.DOES_STACKS;
        ParticleAnalyzer.setRoiManager(countman); 
        int measurements = 0;
        ParticleAnalyzer pa = new ParticleAnalyzer(options, measurements, new ResultsTable(), minSize, Double.MAX_VALUE, 0.0, 1.0);
        
        for (int t=1; t<= nFrames; t++) {
            for (int z=1; z <= nSlices; z++) {
                int index = t*z;
                segmented.setSliceWithoutUpdate(index);
                pa.analyze(segmented);
            }
        }
        
        // Now get the ROIS and rescale them to match the raw image
        Roi[] tmprois = countman.getRoisAsArray();
        countman.runCommand("Deselect");
        countman.runCommand("Delete");

        // FIXME: problem with overlays
        int counter = 1;
        if (display_values) {
            for (int c = firstC; c <= lastC; c++) {
                raw.setC(c); 
                for (Roi tmproi: tmprois) {
                    int pan_x = tmproi.getBounds().x;
                    int pan_y = tmproi.getBounds().y;
                    int new_x = (int) floor(pan_x * invscale);
                    int new_y = (int) floor(pan_y * invscale);
                    
                    Roi tmpscaled = RoiScaler.scale(tmproi, invscale, invscale, false);
                    tmpscaled.setLocation(new_x, new_y);
                    tmpscaled.setName(objName+"-"+counter);
                    raw.setRoi(tmpscaled);
                    measure.measure();
                    singles.addValue("Lysosome Type", objName);
                    singles.addValue("Lysosome Channel", ch_lyso);
                    singles.addValue("Protein Channel", ch_protein);
                    singles.addValue("Measurement Channel", c);
                    singles.show("Results");
                    counter++;
                }
            }
        }

        total = tmprois.length;
        return total;
    }

    /**
     * Same as before, but acts on selected Roi.
     *
     * @see #count(ImagePlus, ImagePlus, int, int, int, String, double, boolean)
     * @param roi to restrict the operations
     * @return int of total count for objClass in roi
     */
    private int count(ImagePlus segmented, ImagePlus raw, Roi roi, int firstC, int lastC, 
                        int objClass, String objName, double minSize, boolean display_values) {
        int width = raw.getWidth();
        int swidth = segmented.getWidth();
        double invscale = (double) width/ (double) swidth;

        int total = 0;

        // Constructors
        RoiManager countman = new RoiManager(true); // Hidden roimanager for this task
        ResultsTable singles = null;
        Frame frame = WindowManager.getFrame("Results");
        if (frame!=null && (frame instanceof TextWindow)) {
            TextWindow tw = (TextWindow)frame;
            ResultsTable table = tw.getTextPanel().getResultsTable();
            if (table!= null) {
                singles = table;
            }
        } else {
                singles = new ResultsTable();
        }
        Analyzer measure = new Analyzer(raw, singles);

        // Get segmented image and apply binary threshold for positives or negatives
        ImageProcessor ip = segmented.getProcessor();
        Calibration cal = segmented.getCalibration();
        ip.setThreshold((double) objClass, (double) objClass, ImageProcessor.NO_LUT_UPDATE);
        segmented.setRoi(roi);
        String roiname = roi.getName();

        // Get the Rois from the segmented image
        double unitSquared = cal.pixelWidth*cal.pixelHeight;
        minSize = minSize / unitSquared; // minsize must be a double in pixel units
        int options = 0;
        options += ParticleAnalyzer.SHOW_NONE;
        options += ParticleAnalyzer.ADD_TO_MANAGER;
        options += ParticleAnalyzer.DOES_STACKS;
        ParticleAnalyzer.setRoiManager(countman); 
        int measurements = 0;
        ParticleAnalyzer pa = new ParticleAnalyzer(options, measurements, new ResultsTable(), minSize, Double.MAX_VALUE, 0.0, 1.0);
        
        for (int t=1; t<= nFrames; t++) {
            for (int z=1; z <= nSlices; z++) {
                int index = t*z;
                segmented.setSliceWithoutUpdate(index);
                pa.analyze(segmented);
            }
        }
       
        // Now get the ROIS and rescale them to match the raw image
        Roi[] tmprois = countman.getRoisAsArray();
        countman.runCommand("Deselect");
        countman.runCommand("Delete");

        // FIXME: problem with overlays
        if (display_values) {
            for (Roi tmproi: tmprois) {
                raw.setZ(tmproi.getZPosition());
                raw.setT(tmproi.getTPosition());
                for (int c = firstC; c <= lastC; c++) {
                    raw.setC(c);

                    int pan_x = tmproi.getBounds().x;
                    int pan_y = tmproi.getBounds().y;
                    int new_x = (int) floor(pan_x * invscale);
                    int new_y = (int) floor(pan_y * invscale);
                    
                    Roi tmpscaled = RoiScaler.scale(tmproi, invscale, invscale, false);
                    tmpscaled.setLocation(new_x, new_y);
                    String tmpname = tmpscaled.getName();
                    tmpscaled.setName(tmpname+"-"+objName);
                    raw.setRoi(tmpscaled);
                    measure.measure();
                    singles.addValue("Lysosome Type", objName);
                    singles.addValue("Lysosome Channel", ch_lyso);
                    singles.addValue("Protein Channel", ch_protein);
                    singles.addValue("Measurement Channel", c);
                    singles.addValue("ROI", roiname);
                    singles.show("Results");
                }
            }
        }

        total = tmprois.length;
        return total;
    }

    /**
     * Pre-processing step. Take multichannel TIFF image and convert it to RGB
     * with Lysosome Channel in green LUT and Protein Channel in red LUT
     *
     * @param ch_protein is the channel of the protein inside lysosomes --> RED
     * @param ch_lyso is the marker for lysosomes --> GREEN
     * @param imp is the TIFF image to convert
     * @return rgb image with settings above
     */
    private ImagePlus make_rgb(int ch_protein, int ch_lyso, ImagePlus imp) {
        Calibration cal = imp.getCalibration();
        ImagePlus[] channels = ChannelSplitter.split(imp);
        IJ.run(channels[ch_protein-1], "Red", "");
        IJ.run(channels[ch_lyso-1], "Green", "");
        
        ImagePlus[] stack = {channels[ch_protein-1], channels[ch_lyso-1]};
        
        ImagePlus mergergb = RGBStackMerge.mergeChannels(stack, false);
        mergergb.setCalibration(cal);
        RGBStackConverter.convertToRGB(mergergb);
        return mergergb;
    }
    
    /**
     * GUI for this plugin
     * 
     * @return true is everything goes fine, false if canceled
     */
    private boolean showDialog() {
	    GenericDialog gd = new GenericDialog("LysoQuant");

	    gd.addNumericField("Lysosome Channel", Integer.parseInt(Prefs.get("lysoquant.display_lyso", "2")), 0);
	    gd.addNumericField("Protein Channel", Integer.parseInt(Prefs.get("lysoquant.display_protein", "3")), 0);

	    if (nSlices > 1) {
            gd.addStringField("Slices", "1-"+nSlices);
	    }	
	
	    if (nFrames > 1) {
            gd.addStringField("Frames", "1-"+nFrames);		      
        }	

	    boolean defaultTick = Boolean.parseBoolean(Prefs.get("lysoquant.display_values", "false"));
	    gd.addCheckbox("Show single values", defaultTick);

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
        
        gd.showDialog();
        if (gd.wasCanceled())
            return false;

        // get entered values
        ch_lyso = (int)gd.getNextNumber();
        ch_protein = (int)gd.getNextNumber();

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
	
        display_values = (boolean)gd.getNextBoolean();
        
        // Retrieve meas. channel range. Taken from Duplicator.java
        String[] range = Tools.split(gd.getNextString(), " -");
        double c1 = gd.parseDouble(range[0]);
        double c2 = range.length==2?gd.parseDouble(range[1]):Double.NaN;
        firstC = Double.isNaN(c1)?1:(int)c1;
        lastC = Double.isNaN(c2)?firstC:(int)c2;
        if (firstC<1) firstC = 1;
        if (lastC>nChannels) lastC = nChannels;
        if (firstC>lastC) {firstC=1; lastC=nChannels;}
        
        // Save for next usage
        Prefs.set("lysoquant.display_lyso", ch_lyso);
        Prefs.set("lysoquant.display_protein", ch_protein);
        Prefs.set("lysoquant.display_values", Boolean.toString(display_values));
        Prefs.set("lysoquant.display_firstC", firstC);
        Prefs.set("lysoquant.display_lastC", lastC);

        if (nSlices > 1 && display_warning) {
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
