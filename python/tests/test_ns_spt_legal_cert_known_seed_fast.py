import json
import subprocess
import sys
from pathlib import Path


def test_ns_spt_legal_cert_known_seed_initial_only(tmp_path):
    root = Path(__file__).resolve().parents[2]
    script = root / "python" / "scripts" / "ns_spt_legal_certify.py"
    out_json = tmp_path / "spt_legal_cert.json"

    proc = subprocess.run(
        [
            sys.executable,
            str(script),
            "--N",
            "16",
            "--seed",
            "0",
            "--t-max",
            "0",
            "--no-evolve",
            "--out-json",
            str(out_json),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    text = (proc.stdout or "") + (proc.stderr or "")
    assert proc.returncode == 0, text
    assert "SPT_LEGAL_VERDICT: PASS" in text
    assert out_json.exists()

    cert = json.loads(out_json.read_text(encoding="utf-8"))
    assert cert["verdict"] == "PASS"
    assert cert["config"]["seed"] == 0
