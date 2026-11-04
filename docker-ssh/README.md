# LysoQuant server in Docker with SSH access

Runs on any Docker host with NVIDIA GPUs: Linux (the usual case), or Windows
with Docker Desktop (WSL2 backend).

The same multi-GPU segmentation server as in [`../server`](../server/README.md),
in a container that also runs its own SSH server (SFTP only). Fiji users on
other computers connect to the container with the **private key** you give
them; the container only trusts the matching **public key**. The SSH accounts
live inside the container, so no account has to be created on the host (for
example on a Windows PC running Docker Desktop).

## Set up (once)

```bash
cd docker-ssh
./setup.sh                 # public key -> authorized_keys; detects NVIDIA or CPU -> .env
docker compose up -d --build
docker compose ps          # "healthy" once the model is loaded
docker compose logs -f
```

`authorized_keys` (one public key per line) is mounted read-only in the
container; `sshd` reads it at every login, so adding or removing a line applies
at once, without a restart (if your editor replaces the file instead of
writing into it, run `docker compose restart`). Use another file with
`LYSOQUANT_AUTHORIZED_KEYS_FILE` in `.env`.

`setup.sh` uses `../key/lysoquant.key` (or a key given as argument) and creates
a new RSA key pair there if the file does not exist; run it again with another
key to add it. The key has to be RSA in
PEM format (`BEGIN RSA PRIVATE KEY`, `ssh-keygen -m PEM`): that is what the SSH
library of Fiji reads. `.env`, `authorized_keys` and `key/` are ignored by git.

`setup.sh` also copies [`env.example`](env.example) to `.env`, where the settings
are, and detects the GPU setup of the computer (`--nvidia` or `--cpu` force one;
run it again after changing the hardware or installing the NVIDIA Container Toolkit):

| Computer | Setup |
|---|---|
| Linux, Windows (Docker Desktop with Git Bash, or WSL) with an NVIDIA GPU that Docker can use | NVIDIA: CUDA image, `LYSOQUANT_GPU=all`, GPU reservation in `compose.override.yaml` |
| No NVIDIA GPU, or Docker without the `nvidia` runtime | CPU: small CPU image, `LYSOQUANT_GPU=none` (slow) |
| macOS | CPU only: Docker on a Mac cannot use the Apple GPU (Metal). For Metal use the native installer in [`../macos`](../macos/README.md) |

`.env` settings: `LYSOQUANT_SSH_PORT` (default **2222**, not 22: it is used
inside the container and on the host, open it in the firewall),
`LYSOQUANT_GPU` (`all` with an NVIDIA GPU, see `../server/README.md`; `none` for
the CPU), `LYSOQUANT_TILE_SIZE` (default 900: size in px of the tiles the images
are segmented in; the largest one that fits in the GPU memory or RAM, with the
estimates in `env.example` and in "Tile size and memory" in `../server/README.md`),
`LYSOQUANT_TILES_PER_GPU` (default 1: tiles segmented at the same time
on each GPU, i.e. jobs run at once per GPU; each needs its own GPU memory, or RAM
on the CPU), `LYSOQUANT_AUTHORIZED_KEYS` (optional
extra public keys, separated by `;`), `LYSOQUANT_JOB_TIMEOUT`,
`LYSOQUANT_KEEP_HOURS`. Requires Docker Compose, and for the NVIDIA setup the
NVIDIA Container Toolkit (on Windows: Docker Desktop with the WSL2 backend and a
current NVIDIA driver). On a Linux host note that Docker-published ports bypass
`ufw` rules (use `iptables` `DOCKER-USER` or bind the port to one address in
`compose.yaml`), and that on SELinux hosts the key file mount needs `:ro,z`.

### Memory (CPU image)

Without a GPU the model runs in RAM: with the default 900 px tile
(`LYSOQUANT_TILE_SIZE`) it needs about 3 GB for the model and one tile, plus the
image. A smaller tile needs less, see the table in `../server/README.md`. Docker Desktop (macOS, and
Windows with the WSL2 backend) caps the memory of all containers: if
`docker inspect <container> --format '{{.State.OOMKilled}}'` prints `true`, or the
log shows `Killed`, raise the limit (Docker Desktop > Settings > Resources, or
`memory=` in `%UserProfile%\.wslconfig`) to at least 6 GB. A job that kills the
server is marked as failed (`error.txt`) after the restart, not retried.

## Fiji (Edit > Options > LysoQuant Settings...)

| Setting | Value |
|---|---|
| Use remote host (SSH) | ticked |
| Server, port | address of the computer running Docker, `2222` (or your port) |
| Username | `lysoquant` |
| RSA key file | the private key, e.g. `key/lysoquant.key` copied to the user's computer |
| Process folder on server | `/jobs` |

"Test connection" in the settings checks the connection, the key and the process
folder, and shows the tile size the server is set to. The tile size is not a Fiji
setting: change `LYSOQUANT_TILE_SIZE` in `.env` and run `docker compose up -d`.

## What is exposed

- sshd accepts only the public keys in `authorized_keys`: no passwords, no root, and only
  SFTP (no shell, no port forwarding). The job folder lives in a Docker volume,
  reachable only through SFTP.
- The segmentation server runs as the unprivileged user `lysoquant`, the same
  account that SFTP logs in as.
- **Everybody with the private key is the same user**: they can see and delete
  each other's job folders (they exist only while a job runs). Whoever gets the
  key file has access, so hand it out only to people who may use the server.
- To change who may connect, edit `authorized_keys`. With one key per person
  (one line each) you can remove a person without changing everybody's key.
- The SSH host keys are kept in a volume, so they survive rebuilds.
