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
import ij.io.FileSaver;
import ij.io.Opener;
import ij.measure.Calibration;

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

    private final String id;
    private final ImagePlus input;

    public Blob(ImagePlus input) {
        this.id = "lysoquant-" + UUID.randomUUID().toString();
        this.input = input;
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
     * Input image as TIFF. Z-slices and frames are flattened, each one is segmented independently.
     */
    public byte[] inputTiff() throws IOException {
        byte[] tiff = new FileSaver(input).serialize();
        if (tiff == null) throw new IOException("Cannot serialize image "+input.getTitle());
        return tiff;
    }

    /**
     * Segmentation parameters as JSON (the tile size is a server setting)
     */
    public byte[] paramsJson() {
        double[] el = getElementSizeUm(input);
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
