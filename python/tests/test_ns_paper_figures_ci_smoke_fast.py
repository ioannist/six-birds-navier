import subprocess
import sys
from pathlib import Path

import numpy as np


def _assert_has_finite_numeric(npz_path: Path) -> None:
    with np.load(npz_path, allow_pickle=False) as data:
        assert "meta_json" in data.files
        found_numeric = False
        for key in data.files:
            arr = data[key]
            if arr.dtype.kind in {"i", "u", "f", "c", "b"} and arr.size > 0:
                vals = np.asarray(arr)
                if np.iscomplexobj(vals):
                    ok = np.all(np.isfinite(np.real(vals))) and np.all(np.isfinite(np.imag(vals)))
                else:
                    ok = np.all(np.isfinite(vals))
                assert ok, f"non-finite values in {npz_path}:{key}"
                found_numeric = True
                break
        assert found_numeric, f"no numeric arrays found in {npz_path}"


def test_ns_paper_figures_ci_smoke_fast(tmp_path):
    root = Path(__file__).resolve().parents[2]
    script = root / "python" / "scripts" / "ns_paper_figures.py"
    proc = subprocess.run(
        [
            sys.executable,
            str(script),
            "--clean",
            "--out-root",
            str(tmp_path),
            "--ids",
            "fig_nsbu12_toy_zeno_gallery",
            "fig_ns_spt_illegal_vs_legal",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    text = (proc.stdout or "") + (proc.stderr or "")
    assert proc.returncode == 0, text
    assert "ns_paper_figures: OK" in text

    png1 = tmp_path / "paper" / "nsbu12_toy_zeno_gallery.png"
    npz1 = tmp_path / "paper" / "nsbu12_toy_zeno_gallery_data.npz"
    png2 = tmp_path / "paper" / "ns_spt_illegal_vs_legal.png"
    npz2 = tmp_path / "paper" / "ns_spt_illegal_vs_legal_data.npz"
    assert png1.exists()
    assert npz1.exists()
    assert png2.exists()
    assert npz2.exists()

    _assert_has_finite_numeric(npz1)
    _assert_has_finite_numeric(npz2)
