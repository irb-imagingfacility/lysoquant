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

import ij.CompositeImage;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.gui.ShapeRoi;
import ij.gui.Toolbar;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.ChannelSplitter;
import ij.plugin.Duplicator;
import ij.plugin.RGBStackConverter;
import ij.plugin.RGBStackMerge;
import ij.plugin.RoiScaler;
import ij.plugin.filter.Analyzer;
import ij.plugin.filter.ParticleAnalyzer;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import ij.process.LUT;

import java.awt.Color;
import java.awt.Rectangle;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import static java.lang.Math.floor;

/**
 * LysoQuant without a user interface: segment the lysosomes of an image on the
 * server, then count them (and measure them, one by one) into result tables.
 * Used by the menu command, the macros, the SciJava commands and scripts.
 * <p>
 * Set the parameters (public fields), then call {@link #analyze} (or
 * {@link #segment} and {@link #count} separately). The tables hold one row per
 * image, slice and frame in {@link #totals}, and one row per lysosome and measurement
 * channel in {@link #singles} (if {@link #singleValues}). Nothing is shown on screen.
 * <p>
 * Channels, slices and frames are counted from 1; a last channel, slice or
 * frame of 0 means up to the end.
 *
 * @author Diego Morone
 */
public class LysoQuantAnalysis {

    /** Class names of the labels in the segmentation, from 1 */
    private static final Map<Integer, String> CLASSES = new LinkedHashMap<>();
    static {
        CLASSES.put(1, "Empty");
        CLASSES.put(2, "Loaded");
    }

    public int lysoChannel = 2;
    public int proteinChannel = 3;
    public int firstC = 1;      // channels to measure
    public int lastC = 0;
    public int firstZ = 1;      // slices to segment
    public int lastZ = 0;
    public int firstT = 1;      // frames to segment
    public int lastT = 0;
    /** Measure every lysosome, in every measurement channel */
    public boolean singleValues = false;
    /** With singleValues: also the % of each lysosome covered by cargo (signal above cargoThreshold) */
    public boolean cargoArea = false;
    public int cargoThreshold = 0;
    /** Smallest lysosome counted, in the squared unit of the image (microns) */
    public double minSize = 0.53;
    /** Draw the measured lysosomes in the overlay of the measured image, as the menu command does */
    public boolean addOverlay = false;

    /** One row per lysosome and measurement channel */
    public final ResultsTable singles;
    /** One row per image, slice and frame: number and fraction of lysosomes of each type */
    public final ResultsTable totals;

    private Overlay overlay;

    public LysoQuantAnalysis() {
        this(new ResultsTable(), new ResultsTable());
    }

    /** Add the results to existing tables */
    public LysoQuantAnalysis(ResultsTable singles, ResultsTable totals) {
        this.singles = singles;
        this.totals = totals;
    }

    /** @return why the image cannot be analyzed, null if it can */
    public static String imageProblem(ImagePlus image) {
        int type = image.getType();
        if (!image.isComposite() && !(type == ImagePlus.GRAY8 || type == ImagePlus.GRAY16)) {
            return "8-bit or 16-bit grayscale stack required";
        }
        if (image.getStackSize() < 2 && !image.isComposite()) {
            return "A minimum of 2 channels are required";
        }
        return null;
    }

    /**
     * Segment the lysosomes of the image on the server, in the slices and frames
     * selected. Cells restrict the segmentation: what is outside of them is erased
     * from the image sent to the server.
     *
     * @param raw image with the lysosome and protein channels
     * @param cells regions of the image to segment, e.g. cells; null for the whole image or its selection
     * @return 16-bit label image (0 background, 1 empty, 2 loaded lysosomes) of the selected slices and frames
     * @throws IllegalArgumentException if the image or the parameters cannot be used
     * @throws IOException if the segmentation fails (message for the user)
     * @throws InterruptedException if canceled
     */
    public ImagePlus segment(ImagePlus raw, Roi[] cells, ServerSettings server) throws IOException, InterruptedException {
        check(raw);
        ImagePlus rgb = makeRgb(raw);
        if (cells != null && cells.length > 0) {
            ShapeRoi union = new ShapeRoi(cells[0]);
            for (int i = 1; i < cells.length; i++) union = union.or(new ShapeRoi(cells[i]));
            clearOutside(rgb, union);
        } else if (raw.getRoi() != null) {
            clearOutside(rgb, raw.getRoi());
        }

        ImagePlus segmented = Segmentation.segment(rgb, server);
        segmented.setDimensions(1, lastZ - firstZ + 1, lastT - firstT + 1);     // only the selected range was segmented
        segmented.setTitle("LQ_" + raw.getTitle());
        return segmented;
    }

    /**
     * Count the lysosomes of a segmentation and add them to the tables
     * ({@link #singles} too, with singleValues).
     *
     * @param segmented result of {@link #segment}
     * @param raw the image that was segmented
     * @param cells the cells segment() was given, to count in each one; null for the whole image
     */
    public void count(ImagePlus segmented, ImagePlus raw, Roi[] cells) {
        check(raw);
        int[] position = {raw.getChannel(), raw.getSlice(), raw.getFrame()};
        Roi selection = raw.getRoi();
        overlay = null;
        try {
            if (cells == null || cells.length == 0) {
                String cellID = "0000-0000";
                if (selection != null) {
                    Rectangle bounds = selection.getBounds();
                    cellID = IJ.pad(bounds.x, 4) + "-" + IJ.pad(bounds.y, 4);
                }
                count(segmented, raw, cellID, null);
            } else {
                double scale = (double) segmented.getWidth() / (double) raw.getWidth();
                for (Roi roi : cells) {
                    int newX = (int) floor(roi.getBounds().x * scale);
                    int newY = (int) floor(roi.getBounds().y * scale);
                    Roi scaled = RoiScaler.scale(roi, scale, scale, false);
                    scaled.setLocation(newX, newY);
                    count(segmented, raw, roi.getName(), scaled);
                }
            }
        } finally {
            if (!addOverlay) {      // leave the image as it was
                raw.setRoi(selection);
                raw.setPosition(position[0], position[1], position[2]);
            }
        }
    }

    /** Segment, then count. See {@link #segment} and {@link #count}. */
    public ImagePlus analyze(ImagePlus raw, Roi[] cells, ServerSettings server) throws IOException, InterruptedException {
        ImagePlus segmented = segment(raw, cells, server);
        count(segmented, raw, cells);
        return segmented;
    }

    /**
     * Check the parameters against the image and complete the ranges (0 or
     * out of range: up to the end of the image).
     *
     * @throws IllegalArgumentException if they cannot be used
     */
    public void check(ImagePlus raw) {
        String problem = imageProblem(raw);
        if (problem != null) throw new IllegalArgumentException(problem);

        int nChannels = raw.getNChannels();
        if (lysoChannel < 1 || lysoChannel > nChannels) {
            throw new IllegalArgumentException("Selected lysosome channel is higher than total number of channels");
        }
        if (proteinChannel < 1 || proteinChannel > nChannels) {
            throw new IllegalArgumentException("Selected protein channel is higher than total number of channels");
        }
        int[] c = range(firstC, lastC, nChannels);
        firstC = c[0]; lastC = c[1];
        int[] z = range(firstZ, lastZ, raw.getNSlices());
        firstZ = z[0]; lastZ = z[1];
        int[] t = range(firstT, lastT, raw.getNFrames());
        firstT = t[0]; lastT = t[1];
    }

    private static int[] range(int first, int last, int n) {
        if (first < 1) first = 1;
        last = (last < 1) ? n : Math.min(last, n);
        return (first > last) ? new int[] {1, n} : new int[] {first, last};
    }

    /**
     * Pre-processing step. Take a multichannel image and convert it to RGB with
     * the lysosome channel in green and the protein channel in red, for the
     * selected slices and frames.
     */
    ImagePlus makeRgb(ImagePlus raw) {
        Calibration cal = raw.getCalibration();
        int nChannels = raw.getNChannels();
        Roi selection = raw.getRoi();
        raw.deleteRoi();         // the whole image is needed
        ImagePlus duplicate;
        try {
            duplicate = new Duplicator().run(raw, 1, nChannels, firstZ, lastZ, firstT, lastT);
        } finally {
            raw.setRoi(selection);
        }
        ImagePlus[] channels = ChannelSplitter.split(duplicate);
        ImagePlus protein = channels[proteinChannel - 1];
        ImagePlus lyso = channels[lysoChannel - 1];
        protein.getProcessor().setColorModel(LUT.createLutFromColor(Color.RED));
        lyso.getProcessor().setColorModel(LUT.createLutFromColor(Color.GREEN));

        ImagePlus merged = RGBStackMerge.mergeChannels(new ImagePlus[] {protein, lyso}, false);
        merged.setCalibration(cal);
        RGBStackConverter.convertToRGB(merged);
        return merged;
    }

    /** Edit>Clear Outside of every slice: the background color outside of the selection */
    private static void clearOutside(ImagePlus image, Roi roi) {
        ImageStack stack = image.getStack();
        for (int i = 1; i <= stack.getSize(); i++) {
            ImageProcessor ip = stack.getProcessor(i);
            ip.setColor(Toolbar.getBackgroundColor());
            ip.fillOutside(roi);
        }
        image.updateImage();
    }

    /**
     * Count the objects of each class in the segmented image, in the selected
     * slices and frames; with singleValues, also measure each object.
     *
     * @param roi optional region (cell) of the segmented image to count in
     */
    private void count(ImagePlus segmented, ImagePlus raw, String cellID, Roi roi) {
        int nSlices = raw.getNSlices();
        int nFrames = raw.getNFrames();
        double invscale = (double) raw.getWidth() / (double) segmented.getWidth();

        if (overlay == null) {
            Overlay existing = raw.getOverlay();
            overlay = (addOverlay && existing != null) ? existing : new Overlay();
        }

        // Particles are collected here instead of in the ROI Manager, which needs a display
        final List<Roi> found = new ArrayList<>();
        Calibration cal = segmented.getCalibration();
        double minPixels = minSize / (cal.pixelWidth * cal.pixelHeight);        // in pixels
        int options = ParticleAnalyzer.SHOW_NONE + ParticleAnalyzer.DOES_STACKS;
        ParticleAnalyzer pa = new ParticleAnalyzer(options, 0, new ResultsTable(), minPixels, Double.MAX_VALUE, 0.0, 1.0) {
            @Override
            protected void saveResults(ImageStatistics stats, Roi particle) {
                found.add((Roi) particle.clone());
            }
        };

        for (int t = firstT; t <= lastT; t++) {
            for (int z = firstZ; z <= lastZ; z++) {
                int[] totalvalues = new int[CLASSES.size()];
                List<Double> cargovalues = new ArrayList<>();

                String slices = (nSlices > 1) ? "-z:" + z + "/" + nSlices : "";
                String frames = (nFrames > 1) ? "-t:" + t + "/" + nFrames : "";
                String roiname = (roi != null && roi.getName() != null) ? roi.getName() : "";
                String title = raw.getTitle() + roiname + slices + frames;

                for (Map.Entry<Integer, String> type : CLASSES.entrySet()) {
                    int objClass = type.getKey();
                    String objName = type.getValue();

                    if (roi != null) segmented.setRoi(roi);

                    // The segmented image starts at firstZ/firstT of the raw image
                    segmented.setT(t - firstT + 1);
                    segmented.setZ(z - firstZ + 1);

                    // Binary threshold for this class
                    ImageProcessor ip = segmented.getProcessor();
                    ip.setThreshold((double) objClass, (double) objClass, ImageProcessor.NO_LUT_UPDATE);
                    found.clear();
                    pa.analyze(segmented);

                    totalvalues[objClass - 1] = found.size();
                    if (!singleValues) continue;

                    // Rescale the objects to match the raw image
                    for (Roi particle : found) {
                        int newX = (int) floor(particle.getBounds().x * invscale);
                        int newY = (int) floor(particle.getBounds().y * invscale);
                        Roi scaled = RoiScaler.scale(particle, invscale, invscale, false);
                        scaled.setLocation(newX, newY);
                        scaled.setName(objName + "-" + (overlay.size() + 1));
                        scaled.setPosition(lysoChannel, z, t);

                        for (int channel = firstC; channel <= lastC; channel++) {
                            raw.setRoi(scaled, false);
                            raw.setPosition(channel, z, t);
                            measure(raw, singles);
                            singles.addValue("Lysosome Type", objName);
                            singles.addValue("Lysosome Channel", lysoChannel);
                            singles.addValue("Protein Channel", proteinChannel);
                            singles.addValue("Measurement Channel", channel);
                            singles.addValue("Image", raw.getTitle());
                            singles.addValue("Cell ID", cellID);
                            if (nSlices > 1) singles.addValue("Slice", slices);
                            if (nFrames > 1) singles.addValue("Frame", frames);

                            // % of the lysosome covered by cargo: signal above cargoThreshold
                            if (cargoArea) {
                                ImageProcessor ipraw = raw.getProcessor();
                                double maxThres = ipraw.getMax();
                                ipraw.setThreshold((double) cargoThreshold, maxThres, ImageProcessor.NO_LUT_UPDATE);
                                // this considers implicitly the ImageStatistics.LIMIT
                                ImageStatistics cargostats = ImageStatistics.getStatistics(ipraw, ImageStatistics.AREA_FRACTION, cal);
                                singles.addValue("%Cargo Area", cargostats.areaFraction);
                                cargovalues.add(cargostats.areaFraction);
                                singles.addValue("Cargo Area minT", (double) cargoThreshold);
                                singles.addValue("Cargo Area maxT", maxThres);
                            }
                        }
                        overlay.add(scaled);
                    }
                }

                updateSummary(title, totalvalues, cargoArea && singleValues ? median(cargovalues) : null);
            }
        }

        if (singleValues && addOverlay) {
            overlay.drawLabels(true);
            overlay.drawNames(true);
            overlay.setLabelColor(Color.white);
            overlay.drawBackgrounds(true);
            raw.setOverlay(overlay);
        }
    }

    /**
     * Measure the selection of the image (what is set in Analyze>Set Measurements...)
     * into a row of the table. A table of its own is used, because ImageJ clears
     * the table it measures into when the Results window is not open.
     */
    private static void measure(ImagePlus image, ResultsTable into) {
        ResultsTable row = new ResultsTable();
        new Analyzer(image, row).measure();
        into.incrementCounter();
        String label = row.getLabel(0);
        if (label != null) into.addLabel(label);
        for (int column = 0; column <= row.getLastColumn(); column++) {
            if (row.columnExists(column)) {
                into.addValue(row.getColumnHeading(column), row.getValueAsDouble(column, 0));
            }
        }
    }

    private static Double median(List<Double> values) {
        if (values.isEmpty()) return Double.NaN;
        Collections.sort(values);
        int n = values.size();
        return (n % 2 == 1) ? values.get(n / 2) : (values.get(n / 2 - 1) + values.get(n / 2)) / 2.0;
    }

    /**
     * Add a row to the summary table
     *
     * @param title of measurement, with image name, roiname, slice and frame
     * @param totalvalues number of objects of each class
     * @param cargomedian median of the % cargo area of the lysosomes, null if not measured
     */
    private void updateSummary(String title, int[] totalvalues, Double cargomedian) {
        int sum = 0;
        for (int n : totalvalues) sum += n;

        totals.incrementCounter();
        totals.addLabel(title);
        totals.addValue("Lysosome Ch", lysoChannel);
        totals.addValue("Protein Ch", proteinChannel);
        for (Map.Entry<Integer, String> type : CLASSES.entrySet()) {
            int n = totalvalues[type.getKey() - 1];
            totals.addValue(type.getValue(), n);
            totals.addValue(type.getValue() + " Ratio", (double) n / (double) sum);
        }
        totals.addValue("Total #", sum);
        if (cargomedian != null) totals.addValue("%Cargo Area Median", cargomedian);
    }
}
