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

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;

/**
 * Minimal file interface to the segmentation server.
 * Implemented over SSH (remote server) or the local file system.
 *
 * @author Diego Morone
 */
public interface Connection {

    /**
     * Write data to path, creating parent folders if needed. The file appears
     * atomically (written as path.part, then renamed) so that a watching server
     * never reads it half-written.
     */
    void put(byte[] data, String path) throws IOException;

    /** Read the whole file at path */
    byte[] get(String path) throws IOException;

    /** True if the file at path exists */
    boolean exists(String path) throws IOException;

    /** Create a folder readable only by its owner (and the server running as root) */
    void createPrivateFolder(String path) throws IOException;

    /** Recursively delete a folder */
    void deleteFolder(String path) throws IOException;

    void disconnect();

    /** Read input stream as a byte array (InputStream.readAllBytes needs Java 9) */
    static byte[] readAllBytes(InputStream inputStream) throws IOException {
        try (InputStream in = inputStream; ByteArrayOutputStream out = new ByteArrayOutputStream()) {
            byte[] buf = new byte[64 * 1024];
            int len;
            while ((len = in.read(buf)) != -1) out.write(buf, 0, len);
            return out.toByteArray();
        }
    }

}
