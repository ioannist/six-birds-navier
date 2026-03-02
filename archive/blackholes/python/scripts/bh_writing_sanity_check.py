#!/usr/bin/env python3
from __future__ import annotations

import re
import sys
from pathlib import Path

INPUT_RE = re.compile(r"\\input\{([^}]+)\}")
LABEL_RE = re.compile(r"\\label\{([^}]+)\}")
REF_RE = re.compile(r"\\(?:ref|eqref|autoref)\{([^}]+)\}")


def resolve_input(base: Path, root: Path, target: str) -> Path:
    p = Path(target)
    if p.suffix == "":
        p = p.with_suffix(".tex")
    if p.is_absolute():
        return p
    candidates = [(base / p).resolve(), (root / p).resolve()]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def scan_tex(root: Path, start: Path) -> dict:
    to_scan = [start]
    visited: set[Path] = set()
    inputs: list[tuple[Path, Path]] = []
    labels: dict[str, list[Path]] = {}
    refs: dict[str, list[Path]] = {}
    missing_inputs: list[Path] = []

    while to_scan:
        path = to_scan.pop()
        if path in visited:
            continue
        visited.add(path)
        try:
            text = path.read_text(encoding="utf-8")
        except FileNotFoundError:
            missing_inputs.append(path)
            continue

        base = path.parent
        for match in INPUT_RE.findall(text):
            target = resolve_input(base, root, match)
            inputs.append((path, target))
            if target not in visited:
                to_scan.append(target)

        for label in LABEL_RE.findall(text):
            labels.setdefault(label, []).append(path)

        for ref in REF_RE.findall(text):
            refs.setdefault(ref, []).append(path)

    # Report missing inputs that were discovered but not found on disk.
    for _, target in inputs:
        if not target.exists() and target not in missing_inputs:
            missing_inputs.append(target)

    return {
        "visited": visited,
        "inputs": inputs,
        "labels": labels,
        "refs": refs,
        "missing_inputs": missing_inputs,
    }


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    paper = root / "docs" / "paper.tex"
    if not paper.exists():
        print(f"Missing {paper}")
        return 1

    data = scan_tex(root, paper)
    labels = data["labels"]
    refs = data["refs"]
    missing_inputs = data["missing_inputs"]

    # Count BH input lines in paper.tex
    paper_text = paper.read_text(encoding="utf-8")
    bh_inputs = [m for m in INPUT_RE.findall(paper_text) if "blackholes" in m]

    duplicate_labels = {k: v for k, v in labels.items() if len(v) > 1}
    missing_refs = sorted([k for k in refs.keys() if k not in labels])

    ok = True
    expected_author_files = [
        root / "tex" / "blackholes" / "01_interface_eft.tex",
        root / "tex" / "blackholes" / "02_conventions.tex",
        root / "tex" / "blackholes" / "03_accounting_constraints.tex",
        root / "tex" / "blackholes" / "04_echo_realizability.tex",
        root / "tex" / "blackholes" / "05_cross_channel_stacking.tex",
        root / "tex" / "blackholes" / "06_protocol_inference.tex",
    ]
    missing_author = [p for p in expected_author_files if not p.exists()]
    if len(bh_inputs) != 1:
        print(f"BH input count in paper.tex: {len(bh_inputs)} (expected 1)")
        ok = False
    else:
        print(f"BH input in paper.tex: {bh_inputs[0]}")

    if missing_inputs:
        ok = False
        print("Missing \\input targets:")
        for path in sorted(set(missing_inputs)):
            print(f"  {path}")
    else:
        print("Missing \\input targets: none")

    if duplicate_labels:
        ok = False
        print("Duplicate labels:")
        for label, files in sorted(duplicate_labels.items()):
            file_list = ", ".join(str(p) for p in files)
            print(f"  {label}: {file_list}")
    else:
        print("Duplicate labels: none")

    if missing_refs:
        ok = False
        print("Missing refs:")
        for ref in missing_refs:
            print(f"  {ref}")
    else:
        print("Missing refs: none")

    if missing_author:
        ok = False
        print("Missing author include files:")
        for path in missing_author:
            print(f"  {path}")
    else:
        print("Missing author include files: none")

    print(f"Scanned {len(data['visited'])} tex file(s).")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
