#!/usr/bin/env python3
from __future__ import annotations

import os
import pathlib
import subprocess
import sys


def main() -> int:
    root = pathlib.Path(__file__).resolve().parents[2]
    script = root / "python" / "scripts" / "ns12_run_all.py"
    env = os.environ.copy()

    include_slow = env.get("NS_INCLUDE_SLOW", "") == "1"
    include_optional = env.get("NS_INCLUDE_OPTIONAL", "") == "1"

    proc = subprocess.run([sys.executable, str(script)], env=env)
    artifacts_dir = root / "python" / "artifacts"
    status = "OK" if proc.returncode == 0 else "FAIL"
    print(
        "ns_run_all summary: "
        f"status={status} "
        f"slow={include_slow} "
        f"optional={include_optional} "
        f"artifacts={artifacts_dir}"
    )
    return proc.returncode


if __name__ == "__main__":
    raise SystemExit(main())
