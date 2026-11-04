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

import java.io.ByteArrayInputStream;
import java.awt.Component;
import java.awt.GraphicsEnvironment;
import java.io.IOException;
import java.util.Vector;

import javax.swing.JOptionPane;
import javax.swing.JPasswordField;

import com.jcraft.jsch.*;

/**
 * SSH/SFTP connection to the segmentation server.
 * Adapted from the Wid-U plugin (Antonello et al., 2023).
 *
 * @author Diego Morone
 */
public class SSHConnection implements Connection {

    private static final int CONNECT_TIMEOUT_MS = 15000;

    private final Session session;

    /**
     * Open a new SSH connection. If keyfile is empty, the user is asked for a password,
     * otherwise the key is used (asking for its passphrase if needed).
     *
     * @param hostname  URL or IP to connect to
     * @param port      SSH port
     * @param username  SSH username
     * @param keyfile   path to private key, or empty for password authentication
     */
    public SSHConnection(String hostname, int port, String username, String keyfile) throws JSchException {
        this(hostname, port, username, keyfile, null);
    }

    /**
     * Same, with the window the password or passphrase prompt belongs to (so that
     * it can be answered while that window is modal).
     */
    public SSHConnection(String hostname, int port, String username, String keyfile, Component parent)
            throws JSchException {
        IJ.showStatus("LysoQuant: connecting to "+hostname+", port "+port+" as "+username);

        JSch jsch = new JSch();
        JSch.setConfig("StrictHostKeyChecking", "no");

        session = jsch.getSession(username, hostname, port);

        // Establish file transfer with compression if possible
        session.setConfig("compression.s2c", "zlib@openssh.com,zlib,none");
        session.setConfig("compression.c2s", "zlib@openssh.com,zlib,none");
        session.setConfig("compression_level", "9");

        if (keyfile == null || keyfile.isEmpty()) {
            session.setPassword(promptPassphrase(parent, "Password for "+username+"@"+hostname));
        } else {
            if (!new java.io.File(keyfile).isFile()) throw new JSchException("SSH key file not found: " + keyfile);
            try {
                jsch.addIdentity(keyfile);
            } catch (JSchException e) {
                // Encrypted key
                jsch.addIdentity(keyfile, promptPassphrase(parent, "Passphrase for key "+keyfile));
            }
        }

        try {
            session.connect(CONNECT_TIMEOUT_MS);
        } catch (JSchException e) {
            session.disconnect();       // a single attempt: do not leave it running in the background
            throw e;
        }
        IJ.showStatus("LysoQuant: SSH connection established");
    }

    @Override
    public void put(byte[] data, String path) throws IOException {
        ChannelSftp sftp = openSftp();
        try {
            mkdirs(sftp, parent(path));
            sftp.put(new ByteArrayInputStream(data), path + ".part");
            sftp.rename(path + ".part", path);
        } catch (SftpException e) {
            throw new IOException("Cannot write "+path+": "+e.getMessage(), e);
        } finally {
            sftp.disconnect();
        }
    }

    @Override
    public byte[] get(String path) throws IOException {
        ChannelSftp sftp = openSftp();
        try {
            return Connection.readAllBytes(sftp.get(path));
        } catch (SftpException e) {
            throw new IOException("Cannot read "+path+": "+e.getMessage(), e);
        } finally {
            sftp.disconnect();
        }
    }

    @Override
    public boolean exists(String path) throws IOException {
        ChannelSftp sftp = openSftp();
        try {
            sftp.lstat(path);
            return true;
        } catch (SftpException e) {
            if (e.id == ChannelSftp.SSH_FX_NO_SUCH_FILE) return false;
            throw new IOException(e.getMessage(), e);
        } finally {
            sftp.disconnect();
        }
    }

    @Override
    public void createPrivateFolder(String path) throws IOException {
        ChannelSftp sftp = openSftp();
        try {
            mkdirs(sftp, parent(path));
            sftp.mkdir(path);
            sftp.chmod(0700, path);
        } catch (SftpException e) {
            throw new IOException("Cannot create "+path+": "+e.getMessage(), e);
        } finally {
            sftp.disconnect();
        }
    }

    @Override
    public void deleteFolder(String path) throws IOException {
        ChannelSftp sftp = openSftp();
        try {
            recursiveFolderDelete(sftp, path);
        } catch (SftpException e) {
            throw new IOException("Cannot delete "+path+": "+e.getMessage(), e);
        } finally {
            sftp.disconnect();
        }
    }

    @Override
    public void disconnect() {
        session.disconnect();
    }

    private ChannelSftp openSftp() throws IOException {
        try {
            ChannelSftp sftp = (ChannelSftp) session.openChannel("sftp");
            sftp.connect(CONNECT_TIMEOUT_MS);
            return sftp;
        } catch (JSchException e) {
            throw new IOException(e.getMessage(), e);
        }
    }

    private static String parent(String path) {
        int i = path.lastIndexOf('/');
        return (i > 0) ? path.substring(0, i) : "/";
    }

    /*
     * Create folder and all missing parents
     */
    private static void mkdirs(ChannelSftp sftp, String path) throws SftpException {
        try {
            sftp.lstat(path);
        } catch (SftpException e) {
            if (e.id != ChannelSftp.SSH_FX_NO_SUCH_FILE) throw e;
            mkdirs(sftp, parent(path));
            sftp.mkdir(path);
        }
    }

    /*
     * Helper function for recursive delete all contents of a remote path
     */
    @SuppressWarnings("unchecked")
    private static void recursiveFolderDelete(ChannelSftp sftp, String path) throws SftpException {
        Vector<ChannelSftp.LsEntry> entries = sftp.ls(path);
        for (ChannelSftp.LsEntry item : entries) {
            String name = item.getFilename();
            if (".".equals(name) || "..".equals(name)) continue;
            if (item.getAttrs().isDir()) {
                recursiveFolderDelete(sftp, path + "/" + name);
            } else {
                sftp.rm(path + "/" + name);
            }
        }
        sftp.rmdir(path);
    }

    /*
     * Ask for password or key passphrase
     */
    private static String promptPassphrase(Component parent, String message) throws JSchException {
        if (GraphicsEnvironment.isHeadless()) {
            throw new JSchException("Cannot ask for a password without a display: set the SSH key file"
                                    + " (without passphrase) in the LysoQuant settings");
        }
        JPasswordField field = new JPasswordField(20);
        int result = JOptionPane.showConfirmDialog(parent, new Object[]{field}, message, JOptionPane.OK_CANCEL_OPTION);
        if (result != JOptionPane.OK_OPTION) throw new JSchException("Authentication canceled");
        return new String(field.getPassword());
    }
}
