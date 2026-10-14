import subprocess
import sys
from pathlib import Path

import pytest


@pytest.mark.slow
def test_ns_spt_legal_sweep_smoke(tmp_path):
    root = Path(__file__).resolve().parents[2]
    script = root / "python" / "scripts" / "ns_spt_legal_sweep.py"
    proc = subprocess.run(
        [
            sys.executable,
            str(script),
            "--resolutions",
            "16",
            "--nu-list",
            "0.002",
            "--seeds",
            "0",
            "--dt",
            "0.001",
            "--t-max",
            "0.05",
            "--record-u-hat-every",
            "5",
            "--out-dir",
            str(tmp_path),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, (proc.stdout or "") + (proc.stderr or "")
    assert (tmp_path / "ns_spt_legal_sweep.png").exists()
    assert (tmp_path / "ns_spt_legal_sweep_data.npz").exists()
