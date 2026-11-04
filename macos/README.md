# LysoQuant server for macOS

Runs the PyTorch segmentation server on the Apple GPU (MPS) as a background
service of your own user: it starts when you log in, and launchd restarts it if
it exits. Everything lives in your home folder, no administrator rights, Homebrew
or Xcode tools are needed. Internet access is needed once (a standalone Python
with [uv](https://github.com/astral-sh/uv) and PyTorch, a few hundred MB).

```bash
./install.sh            # --jobs DIR, --gpu none for the CPU
tail -f ~/Library/Logs/LysoQuant.log     # "Watching ... for jobs" once the model is loaded
```

| What | Where |
|---|---|
| Python, server, model | `~/Library/LysoQuant` |
| Process folder | `~/LysoQuant/jobs` (private, mode 700) |
| Service | `~/Library/LaunchAgents/ch.irb.lysoquant.server.plist` |
| Log | `~/Library/Logs/LysoQuant.log` |

In Fiji (Edit > Options > LysoQuant Settings...) untick "Use remote host" and
set "Process folder on server" to the process folder above (use the full path,
e.g. `/Users/you/LysoQuant/jobs`).

Run `install.sh` again to update or change options. `./uninstall.sh` removes the
service and the installation (`--purge` also deletes the process folder).
launchd only restarts the server when it exits, not when it hangs.

## Making the package

`./build-package.sh` (from the repository) writes `dist/lysoquant-server-macos.zip`
with this installer, the server sources and `model/lysoquant-7-16.pt`. Unzip it
anywhere on the target Mac and run `./install.sh`. The installer also runs
straight from the repository.
