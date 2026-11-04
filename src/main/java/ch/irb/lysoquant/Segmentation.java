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

import ij.IJ;
import ij.ImagePlus;

import java.awt.Component;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.jcraft.jsch.JSchException;

/**
 * Client for the LysoQuant PyTorch segmentation server (server/lysoquant_server/serve.py).
 * <p>
 * Puts the image in a private job folder in the server process folder, waits
 * until the server has segmented it and retrieves the label image. Several
 * users can submit at the same time: the server queues the jobs. Replaces the
 * call to the U-Net Segmentation plugin (caffe_unet backend).
 *
 * @author Diego Morone
 */
public class Segmentation {

    /** Where the server is looked for when the process folder is not set (relative ones from the SSH login directory) */
    private static final String[] USUAL_FOLDERS = {"/jobs", "LysoQuant/jobs"};

    /** Rewritten by the server every second while idle, in the process folder */
    private static final String HEARTBEAT = ".lysoquant-heartbeat.json";

    /** Written by the server in the job folder while segmenting: tiles done and tiles in total */
    private static final String PROGRESS = "progress.txt";

    /** Written by the client when the user cancels: the server stops the job and deletes its folder */
    private static final String CANCEL = "cancel.txt";

    private static final long POLL_MS = 1000;

    /** Esc is looked at more often than the server is polled */
    private static final long ESCAPE_POLL_MS = 100;

    /** How long an idle server may take to rewrite its heartbeat */
    private static final long HEARTBEAT_WAIT_MS = 5000;

    private static final Pattern HEARTBEAT_TIME = Pattern.compile("\"time\":\\s*([0-9.eE+-]+)");
    private static final Pattern PROGRESS_TILES = Pattern.compile("(\\d+)\\s+(\\d+)");
    private static final Pattern HEARTBEAT_JOB = Pattern.compile("\"job\":\\s*(null|\"([^\"]*)\")");
    private static final Pattern HEARTBEAT_DEVICE = Pattern.compile("\"device\":\\s*\"([^\"]*)\"");
    private static final Pattern HEARTBEAT_TILE_SIZE = Pattern.compile("\"tile_size\":\\s*(\\d+)");
    private static final Pattern HEARTBEAT_FOLDER = Pattern.compile("\"folder\":\\s*\"((?:[^\"\\\\]|\\\\.)*)\"");

    /**
     * Segment an RGB image (protein in red, lysosomes in green) with the settings
     * in Edit>Options>LysoQuant Settings..., for the user interface: problems are
     * shown in a dialog.
     *
     * @param rgb image to segment
     * @return 16-bit label image (0 background, 1 empty, 2 loaded), or null on error or if canceled
     */
    public static ImagePlus segment(ImagePlus rgb) {
        try {
            return segment(rgb, ServerSettings.fromPrefs());
        } catch (InterruptedException e) {
            IJ.showStatus("LysoQuant: " + e.getMessage());
            return null;
        } catch (Exception e) {
            IJ.error("LysoQuant", e.getMessage());
            return null;
        } finally {
            IJ.showProgress(1.0);       // hide the progress bar
        }
    }

    /**
     * Segment an RGB image (protein in red, lysosomes in green) on a server. Opens
     * no dialog, so scripts and headless runs can use it. Esc in ImageJ, or an
     * interrupt of the calling thread, cancels the segmentation: the server stops
     * and deletes the job.
     *
     * @param rgb image to segment
     * @param server where the server is
     * @return 16-bit label image (0 background, 1 empty, 2 loaded)
     * @throws IOException with a message for the user if the segmentation fails
     * @throws InterruptedException if canceled
     */
    public static ImagePlus segment(ImagePlus rgb, ServerSettings server) throws IOException, InterruptedException {
        String problem = server.problem();
        if (problem != null) throw new IOException(problem);

        Blob blob = new Blob(rgb);
        String folder = withSlash(server.folder);
        String heartbeat = folder + HEARTBEAT;
        String jobfolder = folder + blob.getID();

        Connection connection = null;
        boolean created = false;
        boolean canceled = false;
        try {
            connection = connect(server, null);
            checkServer(connection, heartbeat, server.remote ? server.hostname : "this computer");

            IJ.showStatus("LysoQuant: sending image");
            connection.createPrivateFolder(jobfolder);
            created = true;
            connection.put(blob.paramsJson(), jobfolder + "/params.json");
            connection.put(blob.inputTiff(), jobfolder + "/input.tif");

            // The server picks up the job once input.tif is there
            while (!connection.exists(jobfolder + "/done.txt")) {
                if (connection.exists(jobfolder + "/error.txt")) {
                    String msg = new String(connection.get(jobfolder + "/error.txt"), StandardCharsets.UTF_8);
                    throw new IOException("Segmentation failed on the server:\n" + msg.trim());
                }
                checkCanceled();
                int[] tiles = readProgress(connection, jobfolder + "/" + PROGRESS);
                if (tiles != null) {
                    IJ.showStatus("LysoQuant: segmenting, tile " + tiles[0] + "/" + tiles[1] + "... (Esc to cancel)");
                    IJ.showProgress(tiles[0], tiles[1]);
                } else {
                    boolean running = heartbeatRunning(connection, heartbeat, blob.getID());
                    IJ.showStatus("LysoQuant: " + (running ? "segmenting" : "waiting in the server queue")
                                  + "... (Esc to cancel)");
                }
                for (long waited = 0; waited < POLL_MS; waited += ESCAPE_POLL_MS) {
                    checkCanceled();
                    try {
                        Thread.sleep(ESCAPE_POLL_MS);
                    } catch (InterruptedException e) {
                        throw new InterruptedException("Canceled by user");
                    }
                }
            }

            IJ.showStatus("LysoQuant: receiving segmentation");
            ImagePlus segmented = blob.segmentation(connection.get(jobfolder + "/output.tif"));
            IJ.showStatus("");
            return segmented;

        } catch (InterruptedException e) {
            canceled = true;
            throw e;
        } finally {
            if (connection != null) {
                if (created && canceled) {
                    // The server may write in the folder while it is deleted and keep it alive:
                    // the marker makes it stop at once and delete the folder itself
                    try {
                        connection.put(new byte[0], jobfolder + "/" + CANCEL);
                    } catch (IOException e) {
                        IJ.log("LysoQuant: cannot tell the server to stop: " + e.getMessage());
                    }
                }
                if (created) {
                    try {
                        connection.deleteFolder(jobfolder);
                    } catch (IOException e) {
                        IJ.log("LysoQuant: cannot delete " + jobfolder + ": " + e.getMessage());
                    }
                }
                connection.disconnect();
            }
        }
    }

    /** Esc in ImageJ, or an interrupt of this thread (e.g. a script that is stopped) */
    private static void checkCanceled() throws InterruptedException {
        if (IJ.escapePressed() || Thread.currentThread().isInterrupted()) {
            throw new InterruptedException("Canceled by user");
        }
    }

    /**
     * @param parent window that asks for the password or the key passphrase, or null
     * @throws IOException with a message for the user if the connection fails
     */
    private static Connection connect(ServerSettings server, Component parent) throws IOException {
        if (!server.remote) return new LocalConnection();
        int sshPort;
        try {
            sshPort = Integer.parseInt(server.port.trim());
        } catch (NumberFormatException e) {
            throw new IOException("The SSH port must be a number: " + server.port);
        }
        try {
            return new SSHConnection(server.hostname, sshPort, server.username, server.keyfile, parent);
        } catch (JSchException e) {
            throw new IOException(sshError(e), e);
        }
    }

    /** What the server reports when it is tested */
    public static class ServerInfo {
        /** e.g. "server running on cuda:0, tile size 900 px" */
        public final String description;
        /** Process folder of the server, as the plugin must address it */
        public final String folder;

        ServerInfo(String description, String folder) {
            this.description = description;
            this.folder = folder;
        }
    }

    /**
     * Check the server settings without segmenting anything: connect and make
     * sure the server is running (the "Test connection" button of the settings).
     * Connects once. If the process folder is empty or holds no server, the usual
     * folders are tried; the server reports its own folder in its heartbeat.
     *
     * @param port   SSH port (ignored for a local server)
     * @param folder process folder set by the user, or empty to look for it
     * @param parent window that asks for the password or the key passphrase
     * @throws IOException with a message for the user if the server cannot be used
     */
    public static ServerInfo test(boolean remote, String hostname, String port, String username, String keyfile,
                                  String folder, Component parent) throws IOException {
        if (remote && hostname.isEmpty()) throw new IOException("Enter the server");

        int sshPort = 22;
        if (remote) {
            try {
                sshPort = Integer.parseInt(port.trim());
            } catch (NumberFormatException e) {
                throw new IOException("Invalid server port: " + port);
            }
        }

        Connection connection = null;
        try {
            connection = remote ? new SSHConnection(hostname, sshPort, username, keyfile, parent) : new LocalConnection();
            String host = remote ? hostname : "this computer";

            List<String> candidates = new ArrayList<>();
            if (!folder.isEmpty()) candidates.add(folder);
            for (String usual : USUAL_FOLDERS) candidates.add(usual);
            if (!remote) candidates.add(System.getProperty("user.home") + "/LysoQuant/jobs");

            String found = null;
            for (String candidate : candidates) {
                try {
                    if (connection.exists(withSlash(candidate) + HEARTBEAT)) {
                        found = candidate;
                        break;
                    }
                } catch (IOException e) {
                    // e.g. no permission on that folder: try the next one
                }
            }
            if (found == null) {
                String where = folder.isEmpty() ? "none of the usual process folders" : folder;
                throw new IOException("The LysoQuant server is not running on " + host + " (no heartbeat in "
                                      + where + ", see server/README.md)");
            }

            String heartbeat = checkServer(connection, withSlash(found) + HEARTBEAT, host);
            Matcher reported = HEARTBEAT_FOLDER.matcher(heartbeat);
            String serverFolder = reported.find() ? reported.group(1).replace("\\\\", "\\") : found;
            return new ServerInfo("server running" + describe(heartbeat), serverFolder);
        } catch (JSchException e) {
            throw new IOException(sshError(e), e);
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new IOException("Interrupted", e);
        } finally {
            if (connection != null) connection.disconnect();
            IJ.showStatus("");
        }
    }

    /** Short message for the common SSH failures */
    private static String sshError(JSchException e) {
        String msg = String.valueOf(e.getMessage());
        if (msg.contains("Auth fail") || msg.contains("USERAUTH fail")) {
            return "Authentication failed: check the username and the SSH key";
        }
        if (msg.contains("UnknownHostException")) return "Unknown host: " + msg.replaceFirst(".*: ", "");
        if (msg.contains("Connection refused")) return "Connection refused: check the server and the port";
        if (msg.contains("timeout")) return "No answer from the server: check the address, the port and the firewall";
        return msg;
    }

    /** Device and tile size of the server, from its heartbeat (older servers do not report the tile size) */
    private static String describe(String heartbeat) {
        Matcher device = HEARTBEAT_DEVICE.matcher(heartbeat);
        Matcher tile = HEARTBEAT_TILE_SIZE.matcher(heartbeat);
        return (device.find() ? " on " + device.group(1) : "")
               + (tile.find() ? ", tile size " + tile.group(1) + " px" : "");
    }

    private static String withSlash(String folder) {
        return folder.endsWith("/") ? folder : folder + "/";
    }

    /**
     * Fail early if the server is not running. A server with a job running
     * counts as running (Docker restarts it if a job hangs); an idle one must
     * rewrite its heartbeat within a few seconds. Only heartbeat values are
     * compared, so the clocks need not be in sync.
     *
     * @return the last heartbeat read
     */
    private static String checkServer(Connection connection, String heartbeat, String host)
            throws IOException, InterruptedException {
        String down = "The LysoQuant server is not running on " + host + " (no heartbeat in "
                      + heartbeat + ", see server/README.md)";
        if (!connection.exists(heartbeat)) throw new IOException(down);

        String first = readHeartbeat(connection, heartbeat);
        if (job(first) != null) return first;
        IJ.showStatus("LysoQuant: contacting server");
        for (long waited = 0; waited < HEARTBEAT_WAIT_MS; waited += POLL_MS / 2) {
            Thread.sleep(POLL_MS / 2);
            String latest = readHeartbeat(connection, heartbeat);
            if (!time(latest).equals(time(first))) return latest;
        }
        throw new IOException(down + ": no heartbeat for " + HEARTBEAT_WAIT_MS / 1000 + " s");
    }

    /** Tiles {done, total} reported by the server, or null if it has not started (or is too old to report) */
    private static int[] readProgress(Connection connection, String progress) {
        try {
            if (!connection.exists(progress)) return null;
            Matcher m = PROGRESS_TILES.matcher(new String(connection.get(progress), StandardCharsets.UTF_8));
            if (!m.find()) return null;
            int total = Integer.parseInt(m.group(2));
            return total > 0 ? new int[] {Integer.parseInt(m.group(1)), total} : null;
        } catch (IOException | NumberFormatException e) {
            return null;        // being replaced, try again at the next poll
        }
    }

    /** Whether the server is working on this job (with several GPUs, several jobs run at once) */
    private static boolean heartbeatRunning(Connection connection, String heartbeat, String id) {
        try {
            return readHeartbeat(connection, heartbeat).contains("\"" + id + "\"");
        } catch (IOException e) {
            return false;       // being replaced, try again at the next poll
        }
    }

    private static String readHeartbeat(Connection connection, String heartbeat) throws IOException {
        return new String(connection.get(heartbeat), StandardCharsets.UTF_8);
    }

    private static String time(String heartbeat) {
        Matcher m = HEARTBEAT_TIME.matcher(heartbeat);
        return m.find() ? m.group(1) : "";
    }

    private static String job(String heartbeat) {
        Matcher m = HEARTBEAT_JOB.matcher(heartbeat);
        return m.find() ? m.group(2) : null;
    }
}
