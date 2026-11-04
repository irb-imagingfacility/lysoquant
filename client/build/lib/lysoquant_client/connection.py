"""Where the process folder of the server is: on this computer or behind SSH (SFTP)."""

import errno
import io
import os
import posixpath
import shutil
import stat


class Connection:
    """The few file operations the protocol needs. Paths use '/'."""

    def exists(self, path):
        raise NotImplementedError

    def get(self, path):
        raise NotImplementedError

    def put(self, data, path):
        """Write a file so that the server never sees it half written"""
        raise NotImplementedError

    def create_private_folder(self, path):
        raise NotImplementedError

    def delete_folder(self, path):
        raise NotImplementedError

    def close(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


class LocalConnection(Connection):
    """The server runs on this computer (or shares its process folder with it)"""

    def exists(self, path):
        return os.path.exists(path)

    def get(self, path):
        with open(path, "rb") as f:
            return f.read()

    def put(self, data, path):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        tmp = path + ".part"
        with open(tmp, "wb") as f:
            f.write(data)
        os.replace(tmp, path)

    def create_private_folder(self, path):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        os.mkdir(path, 0o700)

    def delete_folder(self, path):
        shutil.rmtree(path)


class SFTPConnection(Connection):
    """The server is behind SSH: the files are moved with SFTP, as the Fiji plugin does"""

    def __init__(self, host, port=22, username=None, key_file=None, password=None, timeout=15):
        import paramiko         # only needed for a remote server

        self._ssh = paramiko.SSHClient()
        self._ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())     # as the Fiji plugin
        try:
            self._ssh.connect(host, port=int(port), username=username, key_filename=key_file or None,
                              password=password, timeout=timeout, banner_timeout=timeout, auth_timeout=timeout,
                              look_for_keys=not key_file and password is None, allow_agent=not key_file)
        except paramiko.AuthenticationException as e:
            raise ConnectionError(f"Authentication failed on {host}: check the username and the SSH key") from e
        except OSError as e:        # refused, unknown host, no answer
            raise ConnectionError(f"Cannot connect to {host}:{port}: {e}") from e
        self._sftp = self._ssh.open_sftp()

    def exists(self, path):
        try:
            self._sftp.stat(path)
            return True
        except FileNotFoundError:
            return False

    def get(self, path):
        with self._sftp.open(path, "rb") as f:
            return f.read()

    def put(self, data, path):
        self._makedirs(posixpath.dirname(path))
        tmp = path + ".part"
        self._sftp.putfo(io.BytesIO(data), tmp)
        self._sftp.posix_rename(tmp, path)

    def create_private_folder(self, path):
        self._makedirs(posixpath.dirname(path))
        self._sftp.mkdir(path, 0o700)
        self._sftp.chmod(path, 0o700)

    def delete_folder(self, path):
        for entry in self._sftp.listdir_attr(path):
            child = f"{path}/{entry.filename}"
            if stat.S_ISDIR(entry.st_mode):
                self.delete_folder(child)
            else:
                self._sftp.remove(child)
        self._sftp.rmdir(path)

    def _makedirs(self, path):
        if not path or path == "/" or self.exists(path):
            return
        self._makedirs(posixpath.dirname(path))
        try:
            self._sftp.mkdir(path)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

    def close(self):
        self._sftp.close()
        self._ssh.close()
