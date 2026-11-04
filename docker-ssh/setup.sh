#!/bin/bash
# Prepare this folder for this computer:
#   ./setup.sh [--nvidia | --cpu] [private key]
#
# - authorized_keys: the public key of the private key (default
#   ../key/lysoquant.key; a new RSA key pair in PEM format, the format that Fiji
#   reads, is generated if the file does not exist). Run again to add a key.
# - .env, and compose.override.yaml for the GPU: the GPU setup is detected
#     Linux, Windows (Git Bash, WSL) with an NVIDIA GPU that Docker can use
#                                  -> NVIDIA: CUDA image, all GPUs
#     anything else                -> CPU: small CPU image
#   macOS always gets the CPU: Docker cannot use the Apple GPU (Metal); for that
#   run the native installer in ../macos. --nvidia / --cpu override the detection.
# Run it again after changing the GPU hardware or the Docker installation.
set -euo pipefail
cd "$(dirname "$0")"

MODE=auto
KEY=${HOME}/.ssh/lysoquant.key
while [ $# -gt 0 ]; do
    case "$1" in
        --nvidia) MODE=nvidia ;;
        --cpu) MODE=cpu ;;
        -h | --help) sed -n 2,15p "$0"; exit 0 ;;
        *) KEY="$1" ;;
    esac
    shift
done

# ---- key
if [ ! -f "$KEY" ]; then
    mkdir -p "$(dirname "$KEY")"
    ssh-keygen -q -t rsa -b 4096 -m PEM -N '' -C lysoquant -f "$KEY"
    echo "Generated $KEY (give this file to the Fiji users, keep it private)"
fi
PUBLIC="$(ssh-keygen -y -f "$KEY")" || { echo "Cannot read $KEY (passphrase?)" >&2; exit 1; }
touch authorized_keys
grep -qxF "$PUBLIC" authorized_keys || echo "$PUBLIC" >> authorized_keys
echo "Public key of $KEY is in docker-ssh/authorized_keys"

# ---- platform and GPU
case "$(uname -s)" in
    Darwin) OS=macOS ;;
    MINGW* | MSYS* | CYGWIN*) OS=Windows ;;
    Linux) if grep -qi microsoft /proc/version 2> /dev/null; then OS="Windows (WSL)"; else OS=Linux; fi ;;
    *) OS="$(uname -s)" ;;
esac

detect_gpu() {      # prints nvidia or cpu; the reason goes to stderr
    if [ "$OS" = macOS ]; then
        echo "Docker on macOS cannot use the Apple GPU (Metal): using the CPU. For Metal run ../macos/install.sh" >&2
        echo cpu; return
    fi
    if ! command -v nvidia-smi > /dev/null 2>&1 || ! nvidia-smi -L > /dev/null 2>&1; then
        echo "No NVIDIA GPU found (nvidia-smi): using the CPU" >&2
        echo cpu; return
    fi
    local runtimes
    if runtimes="$(docker info --format '{{json .Runtimes}}' 2> /dev/null)"; then
        case "$runtimes" in
            *nvidia*) ;;
            *) echo "NVIDIA GPU found, but Docker has no nvidia runtime (install the NVIDIA Container Toolkit," \
                    "then run ./setup.sh again): using the CPU" >&2
               echo cpu; return ;;
        esac
    else
        echo "NVIDIA GPU found (could not ask Docker whether it can use it)" >&2
    fi
    echo nvidia
}

if [ "$MODE" = auto ]; then MODE="$(detect_gpu)"; fi

# ---- .env
[ -f .env ] || cp env.example .env
set_env() {         # set_env NAME VALUE: replace the NAME= (or #NAME=) line, or append it
    K="$1" V="$2" awk 'BEGIN { k = ENVIRON["K"]; v = ENVIRON["V"] }
        !done && $0 ~ ("^#?" k "=") { print k "=" v; done = 1; next }
        { print }
        END { if (!done) print k "=" v }' .env > .env.new
    mv .env.new .env
}
CURRENT_GPU="$(grep -E '^LYSOQUANT_GPU=' .env | tail -n 1 | cut -d= -f2- || true)"

if [ "$MODE" = nvidia ]; then
    cp compose.nvidia.yaml compose.override.yaml
    case "${CURRENT_GPU:-none}" in none | cpu) set_env LYSOQUANT_GPU all ;; esac     # keep a custom choice
    set_env LYSOQUANT_BASE_IMAGE pytorch/pytorch:2.5.1-cuda12.4-cudnn9-runtime
    set_env LYSOQUANT_TORCH_INDEX ""
    echo "$OS: NVIDIA GPU setup ($(nvidia-smi -L 2> /dev/null | wc -l | tr -d ' ') GPU(s), CUDA image)"
else
    rm -f compose.override.yaml
    set_env LYSOQUANT_GPU none
    set_env LYSOQUANT_BASE_IMAGE python:3.12-slim
    set_env LYSOQUANT_TORCH_INDEX https://download.pytorch.org/whl/cpu
    echo "$OS: CPU setup (slow: segmenting an image takes minutes)"
fi
chmod 600 .env
echo "Now run: docker compose up -d --build"
