from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


def test_appendix_figure_table_contains_all_manifest_ids() -> None:
    root = Path(__file__).resolve().parents[2]
    script = root / "python" / "scripts" / "ns_paper_appendix_figure_table.py"
    manifest = root / "paper" / "figures_manifest.json"
    output = root / "paper" / "sections" / "A2_figure_table_autogen.tex"

    proc = subprocess.run(
        [sys.executable, str(script)],
        cwd=str(root),
        check=False,
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr

    manifest_entries = json.loads(manifest.read_text(encoding="utf-8"))
    generated = output.read_text(encoding="utf-8")

    for entry in manifest_entries:
        fig_id = entry["id"]
        assert fig_id in generated, f"missing figure id in appendix table: {fig_id}"
