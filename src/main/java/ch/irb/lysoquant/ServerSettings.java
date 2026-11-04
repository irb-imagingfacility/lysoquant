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

import ij.Prefs;

/**
 * Where the segmentation server is and how to reach it. The settings of
 * Edit>Options>LysoQuant Settings... are the default; scripts can pass their own.
 *
 * @author Diego Morone
 */
public class ServerSettings {

    /** true: SSH to hostname, false: the server runs on this computer */
    public boolean remote = true;
    public String hostname = "";
    public String port = "22";
    public String username = "";
    /** Private key file; empty to be asked for a password (not possible without a display) */
    public String keyfile = "";
    /** Process folder of the server, as seen from the SSH login (or from this computer) */
    public String folder = "";

    /** The settings saved by Edit>Options>LysoQuant Settings... */
    public static ServerSettings fromPrefs() {
        ServerSettings s = new ServerSettings();
        s.remote = Boolean.parseBoolean(Prefs.get("lysoquant.remotehost", "true"));
        s.hostname = Prefs.get("lysoquant.server", "");
        s.port = Prefs.get("lysoquant.port", "22");
        s.username = Prefs.get("lysoquant.username", "");
        s.keyfile = Prefs.get("lysoquant.rsaKeyFilename", "");
        s.folder = Prefs.get("lysoquant.processfolder", "");
        return s;
    }

    /**
     * The saved settings with the given values on top of them (null or empty
     * values keep the saved ones), for scripts that name only what differs.
     *
     * @param remote null to keep the saved choice
     */
    public static ServerSettings fromPrefs(Boolean remote, String hostname, String port, String username,
                                           String keyfile, String folder) {
        ServerSettings s = fromPrefs();
        if (remote != null) s.remote = remote;
        if (hostname != null && !hostname.trim().isEmpty()) s.hostname = hostname.trim();
        if (port != null && !port.trim().isEmpty()) s.port = port.trim();
        if (username != null && !username.trim().isEmpty()) s.username = username.trim();
        if (keyfile != null && !keyfile.trim().isEmpty()) s.keyfile = keyfile.trim();
        if (folder != null && !folder.trim().isEmpty()) s.folder = folder.trim();
        return s;
    }

    /** @return a message for the user if the settings cannot be used, null if they are complete */
    public String problem() {
        if ((remote && hostname.isEmpty()) || folder.isEmpty()) {
            return "Please configure the server in Edit>Options>LysoQuant Settings..."
                   + " (the process folder is filled in by \"Test connection\")";
        }
        return null;
    }
}
