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

import org.scijava.command.Command;
import org.scijava.plugin.Parameter;

/**
 * What the LysoQuant commands for scripts (Groovy, Jython, pyimagej,
 * Fiji headless) have in common: the image, the channels, the slices and frames, and
 * where the server is. The server settings default to
 * Edit>Options>LysoQuant Settings...: name only what differs.
 * <p>
 * Channels, slices and frames are counted from 1; 0 as last slice or frame means
 * up to the end. A selection on the image restricts the analysis to it.
 *
 * @author Diego Morone
 */
abstract class LysoQuantCommand implements Command {

    @Parameter(label = "Image", description = "Multichannel image with the lysosome and the protein channels")
    ImagePlus image;

    @Parameter(label = "Lysosome channel", min = "1")
    int lysosomeChannel = 2;

    @Parameter(label = "Protein channel", min = "1")
    int proteinChannel = 3;

    @Parameter(label = "First slice", min = "1")
    int firstSlice = 1;

    @Parameter(label = "Last slice (0: all)", min = "0")
    int lastSlice = 0;

    @Parameter(label = "First frame", min = "1")
    int firstFrame = 1;

    @Parameter(label = "Last frame (0: all)", min = "0")
    int lastFrame = 0;

    @Parameter(label = "Use remote host (SSH)", required = false,
               description = "false: the server runs on this computer. Default: the LysoQuant settings")
    Boolean remote;

    @Parameter(label = "Server", required = false, description = "Default: the LysoQuant settings")
    String server;

    @Parameter(label = "Server port", required = false)
    String port;

    @Parameter(label = "Username", required = false)
    String username;

    @Parameter(label = "SSH key file", required = false,
               description = "Required without a display: a password cannot be asked")
    String keyFile;

    @Parameter(label = "Process folder on server", required = false)
    String processFolder;

    /** The analysis with the parameters of the command */
    LysoQuantAnalysis analysis(LysoQuantAnalysis analysis) {
        analysis.lysoChannel = lysosomeChannel;
        analysis.proteinChannel = proteinChannel;
        analysis.firstZ = firstSlice;
        analysis.lastZ = lastSlice;
        analysis.firstT = firstFrame;
        analysis.lastT = lastFrame;
        return analysis;
    }

    ServerSettings serverSettings() {
        return ServerSettings.fromPrefs(remote, server, port, username, keyFile, processFolder);
    }

    /** Make a failure of the segmentation known to the caller of the command */
    static RuntimeException failure(Exception e) {
        if (e instanceof InterruptedException) {
            Thread.currentThread().interrupt();
            return new RuntimeException("LysoQuant canceled", e);
        }
        return new RuntimeException(e.getMessage(), e);
    }
}
