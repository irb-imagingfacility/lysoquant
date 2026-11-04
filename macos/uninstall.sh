#!/bin/bash
# Remove the LysoQuant server installed by install.sh.
#
#   ./uninstall.sh [--purge]     --purge also deletes the process folder
set -euo pipefail

LABEL=ch.irb.lysoquant.server
PREFIX="$HOME/Library/LysoQuant"
PLIST="$HOME/Library/LaunchAgents/$LABEL.plist"

# The process folder is the last argument of the service
JOBS=""
if [ -f "$PLIST" ]; then
    JOBS="$(/usr/libexec/PlistBuddy -c "Print :ProgramArguments" "$PLIST" | sed '/^[{}]/d;s/^ *//' | tail -n 1)"
fi

launchctl bootout "gui/$(id -u)/$LABEL" 2>/dev/null || true
rm -f "$PLIST" "$HOME/Library/Logs/LysoQuant.log"
rm -rf "$PREFIX"

if [ "${1:-}" = --purge ] && [ -n "$JOBS" ]; then
    rm -rf "$JOBS"
fi
echo "LysoQuant server removed."
