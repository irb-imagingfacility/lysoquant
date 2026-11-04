#!/bin/bash
# Install the LysoQuant server for the current user on macOS. It runs as a
# background service that starts when you log in; no sudo, no Homebrew, no
# Xcode tools (a standalone Python is downloaded with uv).
#
#   ./install.sh [--jobs DIR] [--gpu "GPU 0"|none]
#
# Works from the repository (macos/install.sh, uses ../server and ../model) and
# from the package made by build-package.sh (install.sh next to server/ and model/).
set -euo pipefail

LABEL=ch.irb.lysoquant.server
PREFIX="$HOME/Library/LysoQuant"        # no spaces: the venv scripts have shebangs
PLIST="$HOME/Library/LaunchAgents/$LABEL.plist"
LOG="$HOME/Library/Logs/LysoQuant.log"
PYTHON_VERSION=3.12

HERE="$(cd "$(dirname "$0")" && pwd)"
if [ -d "$HERE/server" ]; then ROOT="$HERE"; else ROOT="$(dirname "$HERE")"; fi
SERVER_DIR="$ROOT/server"
MODEL="$ROOT/model/lysoquant-7-16.pt"
JOBS="$HOME/LysoQuant/jobs"
GPU="GPU 0"       # on a Mac this selects the Apple GPU (MPS); "none" for the CPU

while [ $# -gt 0 ]; do
    case "$1" in
        --jobs) JOBS="$2"; shift 2 ;;
        --gpu)  GPU="$2"; shift 2 ;;
        *) echo "Usage: $0 [--jobs DIR] [--gpu \"GPU 0\"|none]" >&2; exit 2 ;;
    esac
done

[ "$(uname)" = Darwin ] || { echo "This installer is for macOS" >&2; exit 1; }
[ "$(id -u)" != 0 ] || { echo "Run as your own user, without sudo" >&2; exit 1; }
[ -f "$MODEL" ] || { echo "Model $MODEL not found" >&2; exit 1; }
[ -d "$SERVER_DIR" ] || { echo "Server sources $SERVER_DIR not found" >&2; exit 1; }

case "$(uname -m)" in
    arm64)  UV_ARCH=aarch64 ;;
    x86_64) UV_ARCH=x86_64 ;;
    *) echo "Unsupported architecture $(uname -m)" >&2; exit 1 ;;
esac

echo "Installing into $PREFIX"
mkdir -p "$PREFIX/bin" "$PREFIX/models" "$(dirname "$LOG")" "$(dirname "$PLIST")"

if [ ! -x "$PREFIX/bin/uv" ]; then
    echo "Downloading uv"
    TMP="$(mktemp -d)"
    trap 'rm -rf "$TMP"' EXIT
    curl -fsSL "https://github.com/astral-sh/uv/releases/latest/download/uv-$UV_ARCH-apple-darwin.tar.gz" \
        | tar -xz -C "$TMP"
    install -m 755 "$TMP"/uv-*/uv "$PREFIX/bin/uv"
fi

export UV_PYTHON_INSTALL_DIR="$PREFIX/python"
export UV_CACHE_DIR="$PREFIX/cache"
echo "Creating the Python environment (downloads PyTorch, a few hundred MB)"
"$PREFIX/bin/uv" venv --python "$PYTHON_VERSION" --clear "$PREFIX/venv"
"$PREFIX/bin/uv" pip install --python "$PREFIX/venv/bin/python" "$SERVER_DIR"

install -m 644 "$MODEL" "$PREFIX/models/lysoquant-7-16.pt"

# Private to this user: the Fiji plugin of the same user is the only client
mkdir -p "$JOBS"
chmod 700 "$JOBS"

cat > "$PLIST" <<PLIST
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
    <key>Label</key>
    <string>$LABEL</string>
    <key>ProgramArguments</key>
    <array>
        <string>$PREFIX/venv/bin/lysoquant-serve</string>
        <string>--weights</string>
        <string>$PREFIX/models/lysoquant-7-16.pt</string>
        <string>--gpu</string>
        <string>$GPU</string>
        <string>$JOBS</string>
    </array>
    <key>EnvironmentVariables</key>
    <dict>
        <key>PYTHONUNBUFFERED</key>
        <string>1</string>
    </dict>
    <key>RunAtLoad</key>
    <true/>
    <key>KeepAlive</key>
    <true/>
    <key>ThrottleInterval</key>
    <integer>10</integer>
    <key>StandardOutPath</key>
    <string>$LOG</string>
    <key>StandardErrorPath</key>
    <string>$LOG</string>
</dict>
</plist>
PLIST

# Reload if it was already installed
DOMAIN="gui/$(id -u)"
launchctl bootout "$DOMAIN/$LABEL" 2>/dev/null || true
launchctl bootstrap "$DOMAIN" "$PLIST"
launchctl kickstart -k "$DOMAIN/$LABEL"

cat <<MSG

LysoQuant server installed and running (starts at every login).
  Process folder  $JOBS   <- set it as "Process folder on server" in Fiji,
                             with "Use remote host" unticked
  Log             tail -f "$LOG"
  Uninstall       "$HERE/uninstall.sh"
MSG
