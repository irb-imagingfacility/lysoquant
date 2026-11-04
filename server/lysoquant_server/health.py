"""
Health check for lysoquant-serve, e.g. as a Docker HEALTHCHECK.

lysoquant-serve rewrites a heartbeat file in the process folder every poll
(also read by the Fiji plugin to check that the server is up), listing the
jobs being run with their start times. The server is healthy if the heartbeat
is recent and no job has been running for more than the job timeout. Does not
import torch so that it stays fast.

Usage:
    lysoquant-health [--heartbeat FILE] [--max-age S] [--job-timeout S]
"""

import argparse
import json
import os
import sys
import time

HEARTBEAT_NAME = ".lysoquant-heartbeat.json"
HEARTBEAT = os.environ.get("LYSOQUANT_HEARTBEAT", f"/jobs/{HEARTBEAT_NAME}")
MAX_AGE = float(os.environ.get("LYSOQUANT_HEARTBEAT_MAX_AGE", "30"))
JOB_TIMEOUT = float(os.environ.get("LYSOQUANT_JOB_TIMEOUT", "3600"))


def write_heartbeat(path, **state):
    tmp = f"{path}.part"
    with open(tmp, "w") as f:
        json.dump({"time": time.time(), **state}, f)
    os.chmod(tmp, 0o644)        # readable by all users
    os.replace(tmp, path)


def check(path, max_age, job_timeout):
    """Return (healthy, message)."""
    try:
        with open(path) as f:
            state = json.load(f)
    except (OSError, ValueError) as e:
        return False, f"no heartbeat: {e}"
    now = time.time()
    age = now - state["time"]
    if age > max_age:
        return False, f"heartbeat {age:.0f} s old (max {max_age:.0f} s)"
    jobs = state.get("jobs") or {}
    if not jobs:
        return True, f"idle on {state.get('device')}, {state.get('done', 0)} jobs done"
    job, started = min(jobs.items(), key=lambda item: item[1])
    if now - started > job_timeout:
        return False, f"job {job} running for {now - started:.0f} s (timeout {job_timeout:.0f} s)"
    return True, f"{len(jobs)} job(s) running on {state.get('device')}, oldest {job} for {now - started:.0f} s"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--heartbeat", default=HEARTBEAT, help="heartbeat file written by lysoquant-serve")
    parser.add_argument("--max-age", type=float, default=MAX_AGE, help="max heartbeat age when idle (s)")
    parser.add_argument("--job-timeout", type=float, default=JOB_TIMEOUT, help="max duration of one job (s)")
    args = parser.parse_args()
    healthy, msg = check(args.heartbeat, args.max_age, args.job_timeout)
    print(msg)
    sys.exit(0 if healthy else 1)


if __name__ == "__main__":
    main()
