import subprocess
import sys
from pathlib import Path

import pytest


@pytest.mark.slow
def test_nsbu21_route_capacity_numbers_quick(tmp_path):
    root = Path(__file__).resolve().parents[3]
    script = root / "python" / "scripts" / "nsbu21_route_capacity_numbers.py"
    result = subprocess.run(
        [
            sys.executable,
            str(script),
            "--quick",
            "--seed",
            "0",
            "--out-dir",
            str(tmp_path),
        ],
        check=False,
    )
    assert result.returncode == 0

    artifact = tmp_path / "nsbu21_route_capacity_numbers.png"
    assert artifact.exists()
