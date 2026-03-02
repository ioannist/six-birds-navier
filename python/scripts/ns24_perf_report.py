#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "python" / "artifacts"


def _run(cmd: list[str], env: dict[str, str]) -> tuple[int, str, float]:
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, env=env)
    dt = time.perf_counter() - t0
    return proc.returncode, (proc.stdout or "") + (proc.stderr or ""), dt


def _parse_collected(output: str) -> int | None:
    m = re.search(r"(\d+)(?:/\d+)?\s+tests?\s+collected", output)
    return int(m.group(1)) if m else None


def _parse_executed(output: str) -> int:
    counts = 0
    for m in re.finditer(r"(\d+)\s+(passed|failed|error|errors|skipped|xfailed|xpassed)", output):
        counts += int(m.group(1))
    return counts


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="NS CI perf report for fast pytest suite.")
    p.add_argument("--threshold-seconds", type=float, default=5.0)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    ART.mkdir(parents=True, exist_ok=True)

    env = dict(os.environ)
    env["PYTHONPATH"] = str(ROOT / "python" / "src")

    collect_cmd = [sys.executable, "-m", "pytest", "python/tests", "-m", "not slow", "--collect-only"]
    run_cmd = [sys.executable, "-m", "pytest", "python/tests", "-m", "not slow"]

    c_code, c_out, c_wall = _run(collect_cmd, env)
    r_code, r_out, r_wall = _run(run_cmd, env)

    collected = _parse_collected(c_out)
    executed = _parse_executed(r_out)
    threshold_ok = r_wall <= args.threshold_seconds

    print("NS-24 CI perf report")
    print(f"pytest_not_slow_wall={r_wall:.3f}s")
    print(f"collected_not_slow={collected if collected is not None else 'unknown'}")
    print(f"executed_not_slow={executed}")
    print(f"threshold_seconds={args.threshold_seconds:.3f}")
    print(f"threshold_ok={threshold_ok}")

    report = {
        "collect_returncode": c_code,
        "run_returncode": r_code,
        "collect_wall": c_wall,
        "pytest_not_slow_wall": r_wall,
        "collected_not_slow": collected,
        "executed_not_slow": executed,
        "threshold_seconds": args.threshold_seconds,
        "threshold_ok": threshold_ok,
    }
    (ART / "ns24_perf_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    if c_code != 0 or r_code != 0:
        return 1
    return 0 if threshold_ok else 2


if __name__ == "__main__":
    raise SystemExit(main())
