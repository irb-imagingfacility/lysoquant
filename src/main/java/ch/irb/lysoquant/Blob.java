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
import ij.ImageStack;
import ij.io.FileSaver;
import ij.io.Opener;
import ij.measure.Calibration;
import ij.process.ColorProcessor;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.Locale;
import java.util.UUID;

/**
 * A segmentation job: the image to segment, its parameters, and a unique id
 * used as folder name in the server cache folder.
 * <p>
 * Server job folder layout:
 * <pre>
 *   lysoquant-UUID/input.tif     image sent to the server
 *   lysoquant-UUID/params.json   segmentation parameters
 *   lysoquant-UUID/output.tif    16-bit label image (written by the server)
 *   lysoquant-UUID/done.txt      written by the server on success
 *   lysoquant-UUID/error.txt     written by the server on failure
 * </pre>
 *
 * @author Diego Morone
 */
public class Blob {

    /** Pixel size of the model in microns: the image is rescaled to it here, so the server does not have to */
    public static final double MODEL_PIXEL_SIZE_UM = 0.025;

    private final String id;
    private final ImagePlus input;
    private final ImagePlus rescaled;

    /**
     * @param input RGB image, rescaled to the model pixel size (what is sent to the server)
     */
    public Blob(ImagePlus input) {
        this.id = "lysoquant-" + UUID.randomUUID().toString();
        this.input = input;
        this.rescaled = rescale(input);
    }

    /**
     * Simple method to return the associated UUID for a blob
     *
     * @return String   UUID
     */
    public String getID() {
        return id;
    }

    /**
     * Input image as TIFF, at the model pixel size. Z-slices and frames are flattened, each one is segmented independently.
     */
    public byte[] inputTiff() throws IOException {
        byte[] tiff = new FileSaver(rescaled).serialize();
        if (tiff == null) throw new IOException("Cannot serialize image "+input.getTitle());
        return tiff;
    }

    /**
     * Segmentation parameters as JSON (the tile size is a server setting). The image
     * sent is already at the model pixel size, so the server has nothing to rescale.
     */
    public byte[] paramsJson() {
        double[] el = getElementSizeUm(rescaled);
        String json = String.format(Locale.ROOT,
            "{\n  \"element_size_um\": [%s, %s]\n}\n",
            Double.toString(el[0]), Double.toString(el[1]));
        return json.getBytes(StandardCharsets.UTF_8);
    }

    /**
     * Convert the server output into a label image with the calibration of the
     * rescaled input, as the U-Net plugin did.
     */
    public ImagePlus segmentation(byte[] tiff) throws IOException {
        ImagePlus out = new Opener().deserialize(tiff);
        if (out == null) throw new IOException("Cannot read segmentation returned by the server");

        Calibration outcal = out.getCalibration();
        Calibration cal = input.getCalibration().copy();
        cal.setUnit("um");
        cal.pixelWidth = outcal.pixelWidth;
        cal.pixelHeight = outcal.pixelHeight;
        cal.pixelDepth = 1;
        out.setCalibration(cal);
        out.setTitle(input.getTitle() + " - data (segmentation)");
        out.resetDisplayRange();
        return out;
    }

    /**
     * Bilinear rescaling of an RGB image to the model pixel size, slice by slice, as the
     * server did (and the U-Net plugin before it): pixel i of the result is at i / scale in
     * the input, the neighbour beyond the last pixel is mirrored. The result is 8-bit RGB.
     * An image already at the model pixel size, or with no calibration, is returned as is.
     */
    static ImagePlus rescale(ImagePlus rgb) {
        double[] el = getElementSizeUm(rgb);
        double scaleY = el[0] / MODEL_PIXEL_SIZE_UM;
        double scaleX = el[1] / MODEL_PIXEL_SIZE_UM;
        if (scaleY == 1 && scaleX == 1) return rgb;
        if (!(scaleY > 0 && scaleX > 0) || Double.isInfinite(scaleY) || Double.isInfinite(scaleX)) return rgb;

        int w = rgb.getWidth(), h = rgb.getHeight();
        int wOut = (int) Math.floor(w * scaleX + 0.5), hOut = (int) Math.floor(h * scaleY + 0.5);
        int[] xl = new int[wOut], xu = new int[wOut], yl = new int[hOut], yu = new int[hOut];
        float[] dx = new float[wOut], dy = new float[hOut];
        interpolation(w, scaleX, xl, xu, dx);
        interpolation(h, scaleY, yl, yu, dy);

        ImageStack stack = rgb.getStack();
        ImageStack out = new ImageStack(wOut, hOut);
        for (int i = 1; i <= stack.getSize(); i++) {
            int[] in = (int[]) stack.getProcessor(i).convertToRGB().getPixels();
            int[] pixels = new int[wOut * hOut];
            for (int y = 0; y < hOut; y++) {
                int rowLow = yl[y] * w, rowUp = yu[y] * w;
                float fy = dy[y];
                for (int x = 0; x < wOut; x++) {
                    float fx = dx[x];
                    int p00 = in[rowLow + xl[x]], p01 = in[rowLow + xu[x]];
                    int p10 = in[rowUp + xl[x]], p11 = in[rowUp + xu[x]];
                    int value = 0;
                    for (int shift = 16; shift >= 0; shift -= 8) {
                        double top = ((p00 >> shift) & 255) * (1 - fx) + ((p01 >> shift) & 255) * fx;
                        double bottom = ((p10 >> shift) & 255) * (1 - fx) + ((p11 >> shift) & 255) * fx;
                        value |= (int) Math.round(top * (1 - fy) + bottom * fy) << shift;
                    }
                    pixels[y * wOut + x] = value;
                }
            }
            out.addSlice(stack.getSliceLabel(i), new ColorProcessor(wOut, hOut, pixels));
        }

        ImagePlus result = new ImagePlus(rgb.getTitle(), out);
        result.setDimensions(rgb.getNChannels(), rgb.getNSlices(), rgb.getNFrames());
        Calibration cal = rgb.getCalibration().copy();
        cal.setUnit("um");
        cal.pixelWidth = MODEL_PIXEL_SIZE_UM;
        cal.pixelHeight = MODEL_PIXEL_SIZE_UM;
        result.setCalibration(cal);
        return result;
    }

    /** Lower and upper input neighbours of each output pixel along one axis, and the weight of the upper one */
    private static void interpolation(int nIn, double scale, int[] lower, int[] upper, float[] weight) {
        for (int i = 0; i < lower.length; i++) {
            double pos = i / scale;
            int low = (int) Math.floor(pos);
            int up = low + 1;
            if (up >= nIn) up = 2 * (nIn - 1) - up;     // mirrored
            lower[i] = Math.min(low, nIn - 1);
            upper[i] = Math.max(up, 0);
            weight[i] = (float) (pos - low);
        }
    }

    /**
     * Pixel size (y, x) in microns. Same unit handling as the U-Net plugin.
     */
    static double[] getElementSizeUm(ImagePlus imp) {
        Calibration cal = imp.getCalibration();
        double factor = 1;
        switch (cal.getUnit()) {
            case "m": case "meter":       factor = 1000000.0; break;
            case "cm": case "centimeter": factor = 10000.0; break;
            case "mm": case "millimeter": factor = 1000.0; break;
            case "nm": case "nanometer":  factor = 0.001; break;
            case "pm": case "pikometer":  factor = 0.000001; break;
        }
        return new double[] { cal.pixelHeight * factor, cal.pixelWidth * factor };
    }
}
