import json
from pathlib import Path


def test_paper_manifest_exists_and_parses():
    root = Path(__file__).resolve().parents[2]
    manifest = root / "paper" / "figures_manifest.json"
    assert manifest.exists()
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    assert isinstance(payload, list)
    assert len(payload) > 0
