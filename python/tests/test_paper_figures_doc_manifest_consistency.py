import json
from pathlib import Path


def _doc_png_set(figures_md: Path) -> set[str]:
    out: set[str] = set()
    for line in figures_md.read_text(encoding="utf-8").splitlines():
        if not line.startswith("| `python/artifacts/paper/"):
            continue
        parts = line.split("|")
        if len(parts) < 2:
            continue
        first = parts[1].strip()
        if first.startswith("`") and first.endswith("`"):
            out.add(first[1:-1])
    return out


def test_paper_figures_doc_manifest_consistency():
    root = Path(__file__).resolve().parents[2]
    manifest_path = root / "paper" / "figures_manifest.json"
    figures_md = root / "docs" / "fluids" / "blowup" / "figures.md"

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest_pngs = {f"python/artifacts/{entry['artifact_png']}" for entry in manifest}
    doc_pngs = _doc_png_set(figures_md)

    assert manifest_pngs == doc_pngs
