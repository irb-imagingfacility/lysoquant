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

package ch.irb.lysosome_counter;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
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
import ij.process.ImageStatistics;
import static java.lang.Math.floor;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * U-Net Segmentation of Positive and Negative lysosomes.
 * Takes COLOR_RGB images.
 *
 * @author Diego Morone
 */
public class Lysosome_Counter implements PlugIn, Measurements {
        int ch_lyso;
        int ch_protein;
        boolean display_values;
        // positives
        int pos = 2;
        // negatives
        int neg = 1;
        String minSize = Prefs.get("lysoquant.minsize", ""); //in pixel squared
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

            Roi roiA = image.getRoi();
            RoiManager roiman = RoiManager.getInstance();

            int width = image.getWidth();
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
                        Logger.getLogger(Lysosome_Counter.class.getName()).log(Level.SEVERE, null, ex);
                    }

                    ImagePlus segmented = IJ.getImage();
                    
                    String positives = title + "-loaded";
                    count(segmented, positives, pos, minSize, display_values);

                    String negatives = title + "-empty";
                    count(segmented, negatives, neg, minSize, display_values);
                    
                    
                } else {
                    roiman.runCommand("Combine");
                    IJ.run("Clear Outside");
                    
                    try {
                        SegmentationJob.processHyperStack(unetp);
                    } catch (InterruptedException ex) {
                        Logger.getLogger(Lysosome_Counter.class.getName()).log(Level.SEVERE, null, ex);
                    }

                    ImagePlus segmented = IJ.getImage();
                    int swidth = segmented.getWidth();
                    
                    Roi[] rois = roiman.getRoisAsArray();
                    
                    
                    for (Roi roi : rois){
                        
                        String roiName = roi.getName();
                        
                        int pan_x = roi.getBounds().x;
                        int pan_y = roi.getBounds().y;
                        double scale = (double) swidth/ (double) width;
                        
                        int new_x = (int) floor(pan_x * scale);
                        int new_y = (int) floor(pan_y * scale);
                        
                        Roi scaled = RoiScaler.scale(roi, scale, scale, false);
                        scaled.setLocation(new_x, new_y);

                        String positives = title +"-"+ roiName + "-loaded";
                        count(segmented, scaled, positives, pos, minSize, display_values);

                        String negatives = title +"-"+ roiName + "-empty";
                        count(segmented, scaled, negatives, neg, minSize, display_values);
                    }
                    
                }
                
            }
        }     
       
        
        private void count(ImagePlus image, Roi roi, String newImage, int objClass, String size, boolean display_values) {
            //IJ.run("Set Measurements...", "display redirect=None decimal=3");
            ImageProcessor ip = image.getProcessor();
            Calibration cal = image.getCalibration();
            ImageProcessor ip2 = ip.duplicate();
            ImagePlus tmp = new ImagePlus(newImage, ip2);
            tmp.setCalibration(cal);
            tmp.setRoi(roi);
            ip2.setThreshold((double) objClass, (double) objClass, ImageProcessor.NO_LUT_UPDATE);
            String options;
            if (display_values) {
                options = "size="+size+"-Infinity display show=Nothing summarize";
            } else {
                options = "size="+size+"-Infinity show=Nothing summarize";
            }
            IJ.run(tmp, "Analyze Particles...", options);
            
        } 
        
        private void count(ImagePlus image, String newImage, int objClass, ImagePlus originalupscaled, int channel, int size, boolean display_values) {
            //IJ.run("Set Measurements...", "display redirect=None decimal=3");
            int measurements = Analyzer.getMeasurements();
            ImageProcessor ip = image.getProcessor();
            Calibration cal = image.getCalibration();
            ImageProcessor ip2 = ip.duplicate();
            ImagePlus tmp = new ImagePlus(newImage, ip2);
            tmp.setCalibration(cal);
            ip2.setThreshold((double) objClass, (double) objClass, ImageProcessor.NO_LUT_UPDATE);
            int options;
            int maxSize = Integer.MAX_VALUE;
            options = ParticleAnalyzer.SHOW_PROGRESS;
            options += ParticleAnalyzer.SHOW_NONE;
            if (display_values) {
                options += ParticleAnalyzer.SHOW_RESULTS;
            }
            ResultsTable rt = new ResultsTable();
            ParticleAnalyzer pa = new ParticleAnalyzer(options, measurements, rt, size, maxSize);
            if (!pa.analyze(ip2))
                return;
            float[] a = rt.getColumn(ResultsTable.AREA);
            if (a==null)
                return;
            rt.addValue("Mean", stats.mean);
            rt.show("Results");
        } 
        
        private ImagePlus make_rgb(int ch_protein, int ch_lyso, ImagePlus imp) {
            Calibration cal = imp.getCalibration();
            ImagePlus[] channels = ChannelSplitter.split(imp);
            IJ.run(channels[ch_protein-1], "Red", "");
            IJ.run(channels[ch_lyso-1], "Green", "");
            
            ImagePlus[] stack = {channels[ch_protein-1], channels[ch_lyso-1]};
            
            ImagePlus mergergb = RGBStackMerge.mergeChannels(stack, true);
            mergergb.setCalibration(cal);
            RGBStackConverter.convertToRGB(mergergb);
            return mergergb;
        }
        
	private boolean showDialog() {
		GenericDialog gd = new GenericDialog("LysoQuant");

		// default value is 0.00, 2 digits right of the decimal point
		gd.addNumericField("Lysosome Channel", Integer.parseInt(Prefs.get("lysoquant.display_lyso", "2")), 0);
                gd.addNumericField("Protein Channel", Integer.parseInt(Prefs.get("lysoquant.display_protein", "3")), 0);
                gd.addCheckbox("Show single values", Boolean.getBoolean(Prefs.get("lysoquant.display_values", "false")));
                
		gd.showDialog();
		if (gd.wasCanceled())
			return false;

		// get entered values
		ch_lyso = (int)gd.getNextNumber();
                ch_protein = (int)gd.getNextNumber();
                display_values = (boolean)gd.getNextBoolean();
                Prefs.set("lysoquant.display_lyso", ch_lyso);
                Prefs.set("lysoquant.display_protein", ch_protein);
                Prefs.set("lysoquant.display_values", Boolean.toString(display_values));

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
		Class<?> clazz = Lysosome_Counter.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring("file:".length(), url.length() - clazz.getName().length() - ".class".length());
		System.setProperty("plugins.dir", pluginsDir);

		// start ImageJ
		new ImageJ();

		// open the Clown sample
		ImagePlus image = IJ.openImage("~/NetBeansProjects/mef.tif");
		image.show();

		// run the plugin
		IJ.runPlugIn(clazz.getName(), "");
	}
}
