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

import ij.*;
import ij.gui.*;
import ij.Prefs;
import ij.plugin.PlugIn;

/**
 * Dialog for configuring settings for LysoQuant plugin.
 *
 * @author Diego Morone
 */
public class LysoQuantSettings implements PlugIn {

    public final String VERSION="1.0.5";

    @Override
    public void run(String arg) {

        if (!showDialog()) return;

    }

    private boolean showDialog() {

        final String[] gpuList = {
            "none", "all available", "GPU 0", "GPU 1", "GPU 2", "GPU 3",
            "GPU 4", "GPU 5", "GPU 6", "GPU 7" };

        GenericDialog gd = new GenericDialog("LysoQuant Settings...");
        String username;
        String key;
        String gpuflag;
        Boolean remotehost;
        String server;
        String port;
        String processfolder;
        String modelfolder;
        String weights;
        String tilesize;
        Boolean display_warning;

        // Retrieve stored prefs
        if (Prefs.get("lysoquant.username", "") != "") {
            username = Prefs.get("lysoquant.username", "");
        } else {
            username = Prefs.get("unet.username", "");
        }

        if (Prefs.get("lysoquant.rsaKeyFilename", "") != "") {
            key = Prefs.get("lysoquant.rsaKeyFilename", "");
        } else {
            key = Prefs.get("unet.rsaKeyFilename", "");
        }

        if (Prefs.get("lysoquant.modelDefinitionFolder", "") != "") {
            modelfolder = Prefs.get("lysoquant.modelDefinitionFolder", "");
        } else {
            modelfolder = Prefs.get("unet.modelDefinitionFolder", "");
        }
    
        String modelid = Prefs.get("unet.modelId", "");

        if (Prefs.get("lysoquant.weights", "") != "") {
            weights = Prefs.get("lysoquant.weights", "");
        } else if (modelid != "") {
                weights = Prefs.get("unet."+modelid+".weightFile", "");
        } else {
            weights = Prefs.get("lysoquant.weights", "");
        }

        if (Prefs.get("lysoquant.tilesize", "") != "") {
            tilesize = Prefs.get("lysoquant.tilesize", "");
        } else if (modelid != "") {
                tilesize = Prefs.get("unet.segmentation."+modelid+".tileShape_0", "");
        } else {
            tilesize = Prefs.get("lysoquant.tilesize", "");
        }

        String model = Prefs.get("lysoquant.model", "lyso7-16.modeldef.h5");

        if (Prefs.get("lysoquant.gpuflag", "") != "") {
            gpuflag = Prefs.get("lysoquant.gpuflag", "");
        } else {
            gpuflag = Prefs.get("unet.gpuId", "GPU 0");
        }

        if (Prefs.get("lysoquant.remotehost", "") != "") {
            remotehost = Boolean.parseBoolean(Prefs.get("lysoquant.remotehost", ""));
        } else {
            remotehost = Boolean.parseBoolean(Prefs.get("unet.useRemoteHost", "false"));
        }
        
        if (Prefs.get("lysoquant.server", "") != "") {
            server = Prefs.get("lysoquant.server", "");
        } else {
            server = Prefs.get("unet.hostname", "");
        }
        
        if (Prefs.get("lysoquant.port", "") != "") {
            port = Prefs.get("lysoquant.port", "");
        } else {
            port = Prefs.get("unet.port", "");
        }

        if (Prefs.get("lysoquant.processfolder", "") != "") {
            processfolder = Prefs.get("lysoquant.processfolder", "");
        } else {
            processfolder = Prefs.get("unet.processfolder", "");
        }

        if (Prefs.get("lysoquant.display_warning", "") != "") {
            display_warning = Boolean.parseBoolean(Prefs.get("lysoquant.display_warning", ""));
        } else {
            display_warning = true;
        }

        String minsize = Prefs.get("lysoquant.minsize", "0.53");

        // Create interface
        gd.addMessage("LysoQuant - v" + VERSION);
        gd.addStringField("U-Net model folder: ", modelfolder, 50);
        gd.addStringField("U-Net model: ", model, 50);
        gd.addStringField("U-Net model weights: ", weights, 50);
        gd.addCheckbox("U-Net use remote host?: ", remotehost);
        gd.addStringField("U-Net server: ", server, 50);
        gd.addStringField("U-Net server port: ", port, 10);
        gd.addStringField("U-Net username: ", username, 50);
        gd.addStringField("U-Net key: ", key, 50);
        gd.addChoice("U-Net GPU: ", gpuList, gpuflag);
        gd.addStringField("U-Net tile size: ", tilesize);
        gd.addStringField("U-Net process folder: ", processfolder, 50);
        gd.addStringField("Filter min size: ", minsize);
        gd.addCheckbox("Display 3D warning", display_warning);

		gd.showDialog();
		if (gd.wasCanceled())
			return false;

        // Get values from interface
        modelfolder = gd.getNextString();
        model = gd.getNextString();
        weights = gd.getNextString();
        remotehost = gd.getNextBoolean();;
        server = gd.getNextString();
        port = gd.getNextString();
        username = gd.getNextString();
        key = gd.getNextString();
        gpuflag = gd.getNextChoice();
        tilesize = gd.getNextString();
        processfolder = gd.getNextString();
        minsize = gd.getNextString();
        display_warning = gd.getNextBoolean();

        // Store new prefs
        Prefs.set("lysoquant.username", username);
        Prefs.set("lysoquant.rsaKeyFilename", key);
        Prefs.set("lysoquant.modelDefinitionFolder", modelfolder);
        Prefs.set("lysoquant.minsize", minsize);
        Prefs.set("lysoquant.model", model);
        Prefs.set("lysoquant.tilesize", tilesize);
        Prefs.set("lysoquant.weights", weights);
        Prefs.set("lysoquant.gpuflag", gpuflag);
        Prefs.set("lysoquant.remotehost", Boolean.toString(remotehost));
        Prefs.set("lysoquant.server", server);
        Prefs.set("lysoquant.port", port);
        Prefs.set("lysoquant.processfolder", processfolder);
        Prefs.set("lysoquant.display_warning", Boolean.toString(display_warning));

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
