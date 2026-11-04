#!/bin/bash
# Build dist/lysoquant-server-windows.zip: installer, server sources and model in
# one folder, ready to unzip and double-click install.cmd on any Windows PC.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(dirname "$HERE")"
MODEL="$ROOT/model/lysoquant-7-16.pt"
NAME=lysoquant-server-windows
[ -f "$MODEL" ] || { echo "Model $MODEL not found" >&2; exit 1; }

OUT="$ROOT/dist"
STAGE="$(mktemp -d)"
trap 'rm -rf "$STAGE"' EXIT
PKG="$STAGE/$NAME"
mkdir -p "$PKG/model" "$PKG/server"

# Scripts with Windows line endings (cmd.exe is picky about them)
for f in install.ps1 uninstall.ps1 README.md; do
    sed 's/\r$//; s/$/\r/' "$HERE/$f" > "$PKG/$f"
done
cp "$HERE/install.cmd" "$HERE/uninstall.cmd" "$HERE/install-all-users.cmd" "$HERE/uninstall-all-users.cmd" "$PKG/"
cp "$ROOT/server/pyproject.toml" "$ROOT/server/README.md" "$PKG/server/"
cp -R "$ROOT/server/lysoquant_server" "$PKG/server/"
find "$PKG" -name __pycache__ -prune -exec rm -rf {} +
cp "$MODEL" "$PKG/model/"

mkdir -p "$OUT"
rm -f "$OUT/$NAME.zip"
(cd "$STAGE" && zip -qr "$OUT/$NAME.zip" "$NAME")
echo "Built $OUT/$NAME.zip"
