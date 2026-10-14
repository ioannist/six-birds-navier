#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "python" / "src"
sys.path.insert(0, str(SRC))

from nswave.deterministic_npz import sha256_file  # noqa: E402

MANIFEST_PATH = ROOT / "paper" / "figures_manifest.json"
DEFAULT_ARTIFACT_ROOT = ROOT / "python" / "artifacts"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate NS paper figure artifacts from manifest.")
    parser.add_argument("--include-slow", action="store_true")
    parser.add_argument("--include-optional", action="store_true")
    parser.add_argument("--clean", action="store_true")
    parser.add_argument("--verify-determinism", action="store_true")
    parser.add_argument(
        "--out-root",
        default=str(DEFAULT_ARTIFACT_ROOT),
        help="Artifact root directory; outputs are written under <out-root>/paper/...",
    )
    parser.add_argument(
        "--ids",
        nargs="*",
        default=None,
        help="Optional list of figure IDs to run (subset of manifest).",
    )
    return parser.parse_args()


def load_manifest() -> list[dict]:
    return json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))


def select_entries(
    entries: list[dict],
    include_slow: bool,
    include_optional: bool,
    ids: set[str] | None = None,
) -> list[dict]:
    allowed = {"fast"}
    if include_slow:
        allowed.add("slow")
    if include_optional:
        allowed.add("optional")
    selected = [e for e in entries if e["class"] in allowed]
    if ids is not None:
        selected = [e for e in selected if e["id"] in ids]
    return selected


def resolve_out_root(path_str: str) -> Path:
    raw = Path(path_str)
    if raw.is_absolute():
        return raw.resolve()
    return (ROOT / raw).resolve()


def file_path_from_rel(rel: str, artifact_root: Path) -> Path:
    return artifact_root / rel


def clean_manifest_outputs(entries: list[dict], artifact_root: Path) -> None:
    for e in entries:
        for key in ("artifact_png", "artifact_npz"):
            p = file_path_from_rel(e[key], artifact_root)
            if p.exists():
                p.unlink()


def run_entry(
    entry: dict,
    verify_determinism: bool,
    log_lines: list[str],
    artifact_root: Path,
) -> tuple[str, str]:
    script = ROOT / entry["script"]
    out_png = file_path_from_rel(entry["artifact_png"], artifact_root)
    out_npz = file_path_from_rel(entry["artifact_npz"], artifact_root)
    out_dir = out_png.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        sys.executable,
        str(script),
        "--seed",
        str(entry["seed"]),
        "--out-dir",
        str(out_dir),
    ]

    def _exec_once() -> tuple[int, str]:
        proc = subprocess.run(cmd, capture_output=True, text=True)
        text = (proc.stdout or "") + (proc.stderr or "")
        return proc.returncode, text

    code, text = _exec_once()
    log_lines.append(f"[{entry['id']}] RUN1 code={code}\n{text}\n")
    if code != 0:
        return "FAIL", f"nonzero exit ({code})"

    for marker in entry["stdout_markers"]:
        if marker not in text:
            return "FAIL", f"missing marker: {marker}"

    if not out_png.exists():
        return "FAIL", f"missing png: {out_png}"
    if not out_npz.exists():
        return "FAIL", f"missing npz: {out_npz}"

    if verify_determinism:
        hash1 = sha256_file(out_npz)
        code2, text2 = _exec_once()
        log_lines.append(f"[{entry['id']}] RUN2 code={code2}\n{text2}\n")
        if code2 != 0:
            return "FAIL", f"determinism rerun nonzero exit ({code2})"
        for marker in entry["stdout_markers"]:
            if marker not in text2:
                return "FAIL", f"determinism rerun missing marker: {marker}"
        if not out_npz.exists():
            return "FAIL", f"determinism rerun missing npz: {out_npz}"
        hash2 = sha256_file(out_npz)
        if hash1 != hash2:
            return "FAIL", f"npz hash mismatch {hash1} != {hash2}"
        print(f"DETERMINISM {entry['id']}: {hash1} == {hash2}")

    return "OK", "ok"


def main() -> int:
    args = parse_args()
    entries_all = load_manifest()
    selected_ids = None if not args.ids else set(args.ids)
    selected = select_entries(entries_all, args.include_slow, args.include_optional, ids=selected_ids)

    known_ids = {e["id"] for e in entries_all}
    if selected_ids is not None:
        unknown = sorted(selected_ids - known_ids)
        if unknown:
            print(f"unknown ids: {unknown}")
            return 1

    artifact_root = resolve_out_root(args.out_root)
    paper_dir = artifact_root / "paper"
    log_path = paper_dir / "ns_paper_figures.log"

    paper_dir.mkdir(parents=True, exist_ok=True)

    if args.clean:
        # Clean only selected outputs so fast runs don't remove slow paper artifacts
        # referenced by the manuscript.
        clean_manifest_outputs(selected, artifact_root)

    statuses: list[tuple[str, str, str, str, str]] = []
    log_lines: list[str] = []
    overall_ok = True

    for e in selected:
        status, detail = run_entry(e, args.verify_determinism, log_lines, artifact_root)
        if status != "OK":
            overall_ok = False
        statuses.append((e["id"], e["class"], status, e["artifact_png"], e["artifact_npz"]))
        if detail != "ok":
            log_lines.append(f"[{e['id']}] detail: {detail}\n")

    log_path.write_text("\n".join(log_lines), encoding="utf-8")

    print("id | class | status | png | npz")
    for row in statuses:
        print(f"{row[0]} | {row[1]} | {row[2]} | {row[3]} | {row[4]}")
    print(f"ns_paper_figures: {'OK' if overall_ok else 'FAIL'}")
    return 0 if overall_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
