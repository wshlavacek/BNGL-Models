#!/usr/bin/env python3
"""Keep the host awake and stop the Borghans overnight screens by a deadline."""

from __future__ import annotations

import argparse
import datetime as dt
import subprocess
import time
from pathlib import Path


HERE = Path(__file__).resolve().parent


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--deadline", required=True)
    parser.add_argument("--slots", nargs="+", type=int, required=True)
    args = parser.parse_args()

    deadline = dt.datetime.fromisoformat(args.deadline).timestamp()
    log = HERE / "overnight_deadline_guard.log"
    with log.open("a", buffering=1) as stream:
        stream.write(
            f"started\t{dt.datetime.now().astimezone().isoformat(timespec='seconds')}"
            f"\tdeadline={args.deadline}\tslots={args.slots}\n"
        )
        while time.time() < deadline:
            time.sleep(min(60, max(0, deadline - time.time())))

        for slot in args.slots:
            subprocess.run(
                ["screen", "-S", f"borghans-overnight-{slot}", "-X", "quit"],
                check=False,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        stream.write(
            f"stopped\t{dt.datetime.now().astimezone().isoformat(timespec='seconds')}\n"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
