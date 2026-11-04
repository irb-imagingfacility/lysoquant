"""
Segmentation daemon: watch the process folder and run every job dropped in
it by the Fiji clients, one at a time, keeping the model loaded on the GPU.

A job folder lysoquant-<uuid>/ is ready once input.tif and params.json exist
(the client writes them atomically) and neither done.txt nor error.txt does.
Several users can submit at the same time: their jobs are queued. With
several GPUs (--gpu all, or --gpu "GPU 0,GPU 1") each GPU holds its own copy
of the model and takes the oldest waiting job as soon as it is free. Each GPU
segments at most --tiles-per-gpu tiles at the same time (one per running job,
the tiles of a job are run one after the other): GPU memory has to hold that
many tiles. --tile-size sets the size of the tiles: it is a server setting, not a
client one, because it is what determines the GPU (or CPU) memory the server needs.

The server needs read/write access to the process folder and to the job folders
in it. The heartbeat file in the process folder lets clients check that the
server is running. Job folders left over by crashed clients are deleted after
--keep-hours.

Usage:
    lysoquant-serve --weights lysoquant-7-16.pt /path/to/process/folder
"""

import argparse
import os
import queue
import shutil
import signal
import sys
import threading
import time
from pathlib import Path

from .health import HEARTBEAT_NAME, write_heartbeat
from .predict import CANCEL_NAME, DEFAULT_TILE_SIZE, log, prepare_models, run_job


STARTED_NAME = "started.txt"    # present while a job runs: still there after a restart means the server was killed


def fail_interrupted_jobs(root):
    """Mark the jobs that were running when the server was killed (e.g. out of memory) as failed.

    Without it the restarted server would take the same job again and be
    killed again, forever.
    """
    for job in root.glob("lysoquant-*/"):
        try:
            if (job / STARTED_NAME).exists() and not (job / "done.txt").exists() and not (job / "error.txt").exists():
                (job / "error.txt").write_text("The server was killed while processing this job, "
                                               "probably out of memory: give it more memory or use a smaller tile size\n")
                log(f"Job {job.name} was interrupted by a server crash, marked as failed")
        except OSError as e:
            log(f"Cannot mark interrupted job {job.name}: {e}")


def pending_jobs(root):
    """Ready jobs, oldest first. Folders may vanish at any time (client cleanup)."""
    jobs = []
    for job in root.glob("lysoquant-*/"):
        try:
            if ((job / "params.json").is_file() and not (job / "done.txt").exists()
                    and not (job / "error.txt").exists() and not (job / CANCEL_NAME).exists()):
                jobs.append(((job / "input.tif").stat().st_mtime, job))
        except FileNotFoundError:
            pass
    return [job for _, job in sorted(jobs)]


def remove_canceled_jobs(root, running):
    """Delete the waiting jobs the client canceled but could not delete (a running job deletes itself)."""
    for job in root.glob("lysoquant-*/"):
        if job.name not in running and (job / CANCEL_NAME).exists():
            shutil.rmtree(job, ignore_errors=True)
            log(f"Job {job.name} canceled by the client while waiting, deleted")


def remove_stale_jobs(root, keep_hours):
    """Delete job folders not modified for keep_hours (crashed or canceled clients)."""
    limit = time.time() - keep_hours * 3600
    for job in root.glob("lysoquant-*/"):
        try:
            if job.stat().st_mtime < limit:
                shutil.rmtree(job)
                log(f"Removed stale job {job.name}")
        except OSError as e:
            log(f"Cannot remove stale job {job.name}: {e}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("folder", nargs="?", default=os.environ.get("LYSOQUANT_JOBS", "/jobs"),
                        help="process folder to watch (default $LYSOQUANT_JOBS or /jobs)")
    parser.add_argument("--weights", default=os.environ.get("LYSOQUANT_WEIGHTS"),
                        help="model weights (.pt) (default $LYSOQUANT_WEIGHTS)")
    parser.add_argument("--gpu", default=os.environ.get("LYSOQUANT_GPU", "GPU 0"),
                        help="'GPU <n>' or 'none' for CPU (default $LYSOQUANT_GPU or 'GPU 0')")
    parser.add_argument("--tiles-per-gpu", type=int, default=int(os.environ.get("LYSOQUANT_TILES_PER_GPU", "1")),
                        help="tiles segmented at the same time on each GPU (default $LYSOQUANT_TILES_PER_GPU or 1)")
    parser.add_argument("--tile-size", type=int, default=int(os.environ.get("LYSOQUANT_TILE_SIZE", DEFAULT_TILE_SIZE)),
                        help=f"output tile size in px, 4 + 128*k (default $LYSOQUANT_TILE_SIZE or {DEFAULT_TILE_SIZE})")
    parser.add_argument("--interval", type=float, default=float(os.environ.get("LYSOQUANT_POLL", "1")),
                        help="polling interval in seconds (default 1)")
    parser.add_argument("--heartbeat", help=f"heartbeat file (default <folder>/{HEARTBEAT_NAME})")
    parser.add_argument("--keep-hours", type=float, default=float(os.environ.get("LYSOQUANT_KEEP_HOURS", "24")),
                        help="delete job folders untouched for this long (default 24)")
    args = parser.parse_args()
    if not args.weights:
        parser.error("--weights or $LYSOQUANT_WEIGHTS is required")
    if args.tiles_per_gpu < 1:
        parser.error("--tiles-per-gpu must be at least 1")
    if args.tile_size < 1:
        parser.error("--tile-size must be positive")

    root = Path(args.folder)
    if not root.is_dir():
        sys.exit(f"Process folder {root} does not exist")
    heartbeat = args.heartbeat or str(root / HEARTBEAT_NAME)

    stop = []
    signal.signal(signal.SIGTERM, lambda *_: stop.append(True))
    signal.signal(signal.SIGINT, lambda *_: stop.append(True))

    fail_interrupted_jobs(root)
    models = prepare_models(args.weights, args.gpu)
    device = ",".join(str(m[2]) for m in models)
    idle = queue.Queue()            # (model, name) slots not working on a job; the copies of a model share its weights
    for slot in range(args.tiles_per_gpu):
        for model in models:        # one slot per GPU before a second one on any
            idle.put((model, str(model[2]) if args.tiles_per_gpu == 1 else f"{model[2]}#{slot + 1}"))
    lock = threading.Lock()
    running = {}                    # job name -> start time
    workers = []
    done = 0
    aborted = set()                 # jobs whose folder could not be written, do not retry

    def work(job, slot):
        nonlocal done
        try:
            (job / STARTED_NAME).touch()
            if run_job(job, slot[0], tile_size=args.tile_size) == 0:
                with lock:
                    done += 1
        except OSError as e:        # e.g. folder deleted by a canceled client
            log(f"Job {job.name} aborted: {e}")
            with lock:
                aborted.add(job.name)
        finally:
            with lock:
                del running[job.name]
            idle.put(slot)

    def beat():
        with lock:
            jobs = dict(running)
        oldest = min(jobs, key=jobs.get) if jobs else None
        try:
            write_heartbeat(heartbeat, device=device, tile_size=args.tile_size, folder=str(root.resolve()),
                            job=oldest, jobs=jobs, done=done)
        except OSError as e:        # e.g. a client has the file open (Windows)
            log(f"Cannot write the heartbeat: {e}")

    last_cleanup = 0
    log(f"Watching {root} for jobs on {device}, {args.tiles_per_gpu} tile(s) of {args.tile_size} px at a time on each")

    while not stop:
        if not os.access(root, os.R_OK | os.W_OK | os.X_OK):
            log(f"ERROR: cannot read/write {root}")
        else:
            with lock:
                remove_canceled_jobs(root, running)
            for job in pending_jobs(root):
                with lock:
                    if job.name in running or job.name in aborted:
                        continue
                try:
                    slot = idle.get_nowait()
                except queue.Empty:
                    break           # all GPU slots busy, the rest waits in the queue
                with lock:
                    running[job.name] = time.time()
                log(f"Job {job.name} on {slot[1]}")
                worker = threading.Thread(target=work, args=(job, slot), name=slot[1])
                worker.start()
                workers = [w for w in workers if w.is_alive()] + [worker]
            beat()
            if time.time() - last_cleanup > 600:
                remove_stale_jobs(root, args.keep_hours)
                last_cleanup = time.time()
        time.sleep(args.interval)
    for worker in workers:          # finish the running jobs before stopping
        worker.join()
    log("Stopped")


if __name__ == "__main__":
    main()
