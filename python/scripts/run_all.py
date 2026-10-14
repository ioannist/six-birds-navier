import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def main() -> int:
    print("Deprecated: use ns_run_all.py")
    ns_run_all = ROOT / "python" / "scripts" / "ns_run_all.py"
    proc = subprocess.run([sys.executable, str(ns_run_all)])
    return proc.returncode


if __name__ == "__main__":
    raise SystemExit(main())
