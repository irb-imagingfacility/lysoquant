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

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.nio.file.attribute.PosixFilePermissions;
import java.util.Comparator;
import java.util.stream.Stream;

/**
 * Use the segmentation server running on this computer (no SSH).
 *
 * @author Diego Morone
 */
public class LocalConnection implements Connection {

    @Override
    public void put(byte[] data, String path) throws IOException {
        Path p = Paths.get(path);
        Files.createDirectories(p.getParent());
        Path part = Paths.get(path + ".part");
        Files.write(part, data);
        Files.move(part, p, StandardCopyOption.REPLACE_EXISTING, StandardCopyOption.ATOMIC_MOVE);
    }

    @Override
    public byte[] get(String path) throws IOException {
        return Files.readAllBytes(Paths.get(path));
    }

    @Override
    public boolean exists(String path) {
        return new File(path).exists();
    }

    @Override
    public void createPrivateFolder(String path) throws IOException {
        Path p = Paths.get(path);
        Files.createDirectories(p.getParent());
        Files.createDirectory(p);
        try {
            Files.setPosixFilePermissions(p, PosixFilePermissions.fromString("rwx------"));
        } catch (UnsupportedOperationException e) {
            // not a POSIX file system (Windows)
        }
    }

    @Override
    public void deleteFolder(String path) throws IOException {
        try (Stream<Path> files = Files.walk(Paths.get(path))) {
            for (Path f : (Iterable<Path>) files.sorted(Comparator.reverseOrder())::iterator) {
                Files.delete(f);
            }
        }
    }

    @Override
    public void disconnect() {
    }
}
