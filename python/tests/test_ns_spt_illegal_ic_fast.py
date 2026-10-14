import re
import subprocess
import sys
from pathlib import Path


def _parse_case_line(text: str, case: str) -> tuple[str, float, float]:
    pat = rf"CASE={case}\s+verdict=(PASS|FAIL)\s+max_kappa=([0-9.eE+-]+)\s+max_burst=([0-9.eE+-]+)"
    m = re.search(pat, text)
    if m is None:
        raise AssertionError(f"missing CASE line for {case}:\n{text}")
    return m.group(1), float(m.group(2)), float(m.group(3))


def test_ns_spt_illegal_ic_demo_fast(tmp_path):
    root = Path(__file__).resolve().parents[2]
    script = root / "python" / "scripts" / "ns_spt_illegal_ic_demo.py"
    proc = subprocess.run(
        [
            sys.executable,
            str(script),
            "--N",
            "16",
            "--seed",
            "0",
            "--j-target",
            "3",
            "--window-c",
            "0.5",
            "--no-evolve",
            "--out-dir",
            str(tmp_path),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    text = (proc.stdout or "") + (proc.stderr or "")
    assert proc.returncode == 0, text

    legal_verdict, legal_kappa, _ = _parse_case_line(text, "LEGAL")
    illegal_verdict, illegal_kappa, _ = _parse_case_line(text, "ILLEGAL")

    assert legal_verdict == "PASS"
    assert illegal_verdict == "FAIL"
    assert legal_kappa < 10.0
    assert illegal_kappa > 50.0

    assert (tmp_path / "ns_spt_illegal_vs_legal.png").exists()
    assert (tmp_path / "ns_spt_illegal_vs_legal_data.npz").exists()
