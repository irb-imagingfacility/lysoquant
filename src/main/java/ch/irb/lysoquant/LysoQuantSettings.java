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

import java.awt.Button;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.Insets;
import java.awt.Panel;
import java.awt.TextField;
import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;
import java.util.Vector;

import javax.swing.JLabel;
import javax.swing.SwingUtilities;

/**
 * Dialog for configuring settings for LysoQuant plugin.
 *
 * @author Diego Morone
 */
public class LysoQuantSettings implements PlugIn {

    /** Version of the pom this was built from (lysoquant.properties is filtered by Maven) */
    public static final String VERSION = readVersion();

    private static String readVersion() {
        try (InputStream in = LysoQuantSettings.class.getResourceAsStream("/lysoquant.properties")) {
            if (in != null) {
                Properties properties = new Properties();
                properties.load(in);
                String version = properties.getProperty("version", "").trim();
                if (!version.isEmpty() && !version.startsWith("${")) return version;
            }
        } catch (IOException e) {
            // shown as unknown
        }
        return "unknown";   // e.g. not built with Maven, which fills in the version
    }

    @Override
    public void run(String arg) {

        if (!showDialog()) return;

    }

    /*
     * LysoQuant preference, falling back to the old U-Net plugin preference
     */
    private static String pref(String key, String unetKey, String def) {
        String value = Prefs.get("lysoquant."+key, "");
        if (!value.isEmpty()) return value;
        return (unetKey != null) ? Prefs.get("unet."+unetKey, def) : def;
    }

    private boolean showDialog() {

        GenericDialog gd = new GenericDialog("LysoQuant Settings...");

        // Retrieve stored prefs
        boolean remotehost = Boolean.parseBoolean(pref("remotehost", "useRemoteHost", "true"));
        String server = pref("server", "hostname", "");
        String port = pref("port", "port", "22");
        String username = pref("username", "username", "");
        String key = pref("rsaKeyFilename", "rsaKeyFilename", "");
        String processfolder = pref("processfolder", null, "");
        String minsize = pref("minsize", null, "0.53");
        boolean display_warning = Boolean.parseBoolean(pref("display_warning", null, "true"));
        boolean hide_segmentation = Boolean.parseBoolean(pref("hide_segmentation", null, "false"));

        // Create interface
        gd.addMessage("LysoQuant - v" + VERSION);
        gd.addCheckbox("Use remote host (SSH)", remotehost);
        gd.addStringField("Server: ", server, 50);
        gd.addStringField("Server port: ", port, 10);
        gd.addStringField("Username: ", username, 50);
        gd.addStringField("SSH key (empty for password): ", key, 50);
        addTestConnection(gd);
        gd.addStringField("Process folder on server: ", processfolder, 50);
        gd.addStringField("Filter min size: ", minsize);
        gd.addCheckbox("Display 3D warning", display_warning);
        gd.addCheckbox("Hide the segmentation image (LQ)", hide_segmentation);

		gd.showDialog();
		if (gd.wasCanceled())
			return false;

        // Get values from interface
        remotehost = gd.getNextBoolean();
        server = gd.getNextString();
        port = gd.getNextString();
        username = gd.getNextString();
        key = gd.getNextString();
        processfolder = gd.getNextString();
        minsize = gd.getNextString();
        display_warning = gd.getNextBoolean();
        hide_segmentation = gd.getNextBoolean();

        // Store new prefs
        Prefs.set("lysoquant.remotehost", Boolean.toString(remotehost));
        Prefs.set("lysoquant.server", server);
        Prefs.set("lysoquant.port", port);
        Prefs.set("lysoquant.username", username);
        Prefs.set("lysoquant.rsaKeyFilename", key);
        Prefs.set("lysoquant.processfolder", processfolder);
        Prefs.set("lysoquant.minsize", minsize);
        Prefs.set("lysoquant.display_warning", Boolean.toString(display_warning));
        Prefs.set("lysoquant.hide_segmentation", Boolean.toString(hide_segmentation));

        return true;
    }

    /**
     * "Test connection" button with its result next to it, to be added after the
     * connection fields. It tests the values currently in the dialog (the fields
     * must be added in the order of showDialog), without storing them.
     */
    private void addTestConnection(final GenericDialog gd) {
        final Button button = new Button("Test connection");
        final JLabel status = new JLabel(" ");
        // The preferred height of the button is 0 before the dialog is shown (the text would be
        // invisible): use that of the label, at least as high as a button to be centered on it
        status.setPreferredSize(new Dimension(380, Math.max(status.getPreferredSize().height, 29)));

        Panel panel = new Panel(new FlowLayout(FlowLayout.LEFT, 8, 0));
        panel.add(button);
        panel.add(status);
        gd.addPanel(panel, GridBagConstraints.WEST, new Insets(0, 0, 5, 0));

        button.addActionListener(e -> {
            // Same order as the fields of the dialog
            Vector<?> text = gd.getStringFields();
            final boolean remote = ((java.awt.Checkbox) gd.getCheckboxes().get(0)).getState();
            final String server = ((TextField) text.get(0)).getText().trim();
            final String port = ((TextField) text.get(1)).getText().trim();
            final String username = ((TextField) text.get(2)).getText().trim();
            final String key = ((TextField) text.get(3)).getText().trim();
            final TextField folderField = (TextField) text.get(4);
            final String folder = folderField.getText().trim();

            button.setEnabled(false);
            showStatus(status, "Testing...", Color.GRAY, null);
            Thread tester = new Thread(() -> {
                try {
                    Segmentation.ServerInfo info = Segmentation.test(remote, server, port, username, key, folder, gd);
                    SwingUtilities.invokeLater(() -> folderField.setText(info.folder));
                    showStatus(status, "Connected: " + info.description, new Color(0, 130, 0), null);
                } catch (Exception ex) {
                    String message = (ex.getMessage() != null) ? ex.getMessage() : ex.toString();
                    showStatus(status, message, new Color(190, 0, 0), message);
                } finally {
                    SwingUtilities.invokeLater(() -> button.setEnabled(true));
                }
            }, "LysoQuant connection test");
            tester.setDaemon(true);
            tester.start();
        });
    }

    /** Show a one-line result (the first line of text, shortened; the full text as tooltip) */
    private static void showStatus(final JLabel label, String text, final Color color, final String tooltip) {
        String line = text.split("\n", 2)[0];
        final String shown = (line.length() > 70) ? line.substring(0, 67) + "..." : line;
        SwingUtilities.invokeLater(() -> {
            label.setText(shown);
            label.setForeground(color);
            label.setToolTipText(tooltip);
        });
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
