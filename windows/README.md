
# LysoQuant server for Windows

Runs the PyTorch segmentation server as a background task of your own user: it
starts when you log in and Task Scheduler restarts it if it crashes. Everything
lives in your user folder; no administrator rights and no Python are needed.
Internet access is needed once (a standalone Python with
[uv](https://github.com/astral-sh/uv) and PyTorch: a few hundred MB, about 2.5 GB
with CUDA).

Unzip the package anywhere and double-click `install.cmd` (or run
`.\install.ps1 [-Jobs DIR] [-Gpu none]` in PowerShell). With an NVIDIA GPU and
driver the CUDA version of PyTorch is installed; otherwise the server runs on
the CPU, which is slow (use `-Gpu none` to force it, `-Gpu all` to use all GPUs).

| What | Where |
|---|---|
| Python, server, model | `%LOCALAPPDATA%\LysoQuant` |
| Process folder | `%USERPROFILE%\LysoQuant\jobs` |
| Task | "LysoQuant Server" in Task Scheduler |
| Log | `%LOCALAPPDATA%\LysoQuant\LysoQuant.log` |

In Fiji (Edit > Options > LysoQuant Settings...) untick "Use remote host" and
set "Process folder on server" to the process folder above, e.g.
`C:/Users/you/LysoQuant/jobs`.

Run the installer again to update or change options. `uninstall.cmd` removes the
task and the installation (`-Purge` also deletes the process folder).

## One server for all users of the PC

Double-click `install-all-users.cmd` (it asks for administrator rights; or run
`.\install.ps1 -AllUsers` in an administrator PowerShell). The server is
installed in `C:\Program Files\LysoQuant`, runs as SYSTEM from startup (no one
needs to be logged in) as the task "LysoQuant Server (all users)", and all
users share one process folder, `C:\ProgramData\LysoQuant\jobs`, and one GPU:
jobs are queued and run one at a time per GPU (`-Gpu all` uses every NVIDIA GPU). Log: `C:\ProgramData\LysoQuant\LysoQuant.log`.
Every user sets in Fiji (LysoQuant Settings, "Use remote host" unticked) the
process folder `C:/ProgramData/LysoQuant/jobs`. Remove it with
`uninstall-all-users.cmd`. Do not install the per-user and the all-users version
on the same PC.

Users can read and delete each other's job folders in this folder (the files
live there only while a job runs).

## Making the package

`./build-package.sh` (from the repository, on macOS or Linux) writes
`dist/lysoquant-server-windows.zip`.
