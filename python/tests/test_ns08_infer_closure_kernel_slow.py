import subprocess
import sys
from pathlib import Path

import pytest


@pytest.mark.slow
def test_ns08_infer_closure_kernel_runs():
    root = Path(__file__).resolve().parents[2]
    script = root / "python" / "scripts" / "ns08_infer_closure_kernel.py"

    proc = subprocess.run(
        [sys.executable, str(script)],
        capture_output=True,
        text=True,
        check=False,
    )

    assert proc.returncode == 0, proc.stdout + proc.stderr

    artifact = root / "python" / "artifacts" / "ns08_closure_kernel.png"
    assert artifact.exists()
    assert "closure_dissipative=True" in proc.stdout
