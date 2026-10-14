#!/usr/bin/env bash
set -euo pipefail

# Create a zip of repo sources plus reproducibility artifacts.
# Artifacts under python/artifacts are always included.
repo_root="$(cd "$(dirname "$0")/.." && pwd)"
REPO_ROOT="$repo_root" python3 - <<'PY'
from __future__ import annotations

from dataclasses import dataclass
import fnmatch
import os
import pathlib
import zipfile

root = pathlib.Path(os.environ["REPO_ROOT"]).resolve()
version_file = root / ".package-repo-snapshot-version"
allowed_roots = [root / "python", root / "formal", root / "docs", root / "scripts", root / "paper"]
allowed_root_files = {"check_python.sh", "check_lean.sh", "check_paper.sh"}
excluded_dirs = {
    ".git",
    ".lake",
    ".venv",
    "venv",
    "__pycache__",
    ".mypy_cache",
    ".pytest_cache",
    ".ruff_cache",
}
allowed_python_ext = {".py", ".txt", ".ini", ".toml", ".png"}
allowed_scripts_ext = {".sh", ".py", ".txt"}
allowed_formal_names = {"lean-toolchain", "lakefile.lean", "lakefile.toml"}

current_version = -1
if version_file.exists():
    raw = version_file.read_text(encoding="utf-8").strip()
    if raw.isdigit():
        current_version = int(raw)

next_version = current_version + 1
zip_path = root / f"repo_snapshot_v{next_version}.zip"

@dataclass(frozen=True)
class IgnoreRule:
    base: pathlib.Path
    pattern: str
    negated: bool
    anchored: bool
    dir_only: bool


def _parse_gitignore(path: pathlib.Path) -> list[IgnoreRule]:
    rules: list[IgnoreRule] = []
    if not path.exists():
        return rules
    base = path.parent
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith(r"\#") or line.startswith(r"\!"):
            line = line[1:]
        negated = line.startswith("!")
        if negated:
            line = line[1:]
        anchored = line.startswith("/")
        if anchored:
            line = line[1:]
        dir_only = line.endswith("/")
        if dir_only:
            line = line[:-1]
        if not line:
            continue
        rules.append(
            IgnoreRule(
                base=base,
                pattern=line,
                negated=negated,
                anchored=anchored,
                dir_only=dir_only,
            )
        )
    return rules


def _match_rule(rule: IgnoreRule, rel_posix: str, is_dir: bool) -> bool:
    if rule.dir_only and not is_dir:
        # Directory rules apply to directories and their contents.
        parts = rel_posix.split("/")
        for i in range(1, len(parts) + 1):
            prefix = "/".join(parts[:i])
            if _match_rule(
                IgnoreRule(
                    base=rule.base,
                    pattern=rule.pattern,
                    negated=rule.negated,
                    anchored=rule.anchored,
                    dir_only=False,
                ),
                prefix,
                is_dir=True,
            ):
                return True
        return False
    if rule.anchored:
        return fnmatch.fnmatchcase(rel_posix, rule.pattern)
    if "/" in rule.pattern:
        return fnmatch.fnmatchcase(rel_posix, rule.pattern)
    return fnmatch.fnmatchcase(pathlib.PurePosixPath(rel_posix).name, rule.pattern)


def _collect_rules(root_path: pathlib.Path) -> list[IgnoreRule]:
    rules: list[IgnoreRule] = []
    rules.append(
        IgnoreRule(
            base=root_path,
            pattern=".git",
            negated=False,
            anchored=False,
            dir_only=True,
        )
    )
    rules.append(
        IgnoreRule(
            base=root_path,
            pattern="repo_snapshot.zip",
            negated=False,
            anchored=False,
            dir_only=False,
        )
    )
    rules.append(
        IgnoreRule(
            base=root_path,
            pattern="repo_snapshot_v*.zip",
            negated=False,
            anchored=False,
            dir_only=False,
        )
    )
    for dirpath, dirnames, filenames in os.walk(root_path):
        if ".git" in dirnames:
            dirnames.remove(".git")
        ignore_path = pathlib.Path(dirpath) / ".gitignore"
        if ignore_path.exists():
            rules.extend(_parse_gitignore(ignore_path))
    return rules


def _is_ignored(path: pathlib.Path, rules: list[IgnoreRule]) -> bool:
    rel_posix = path.relative_to(root).as_posix()
    ignored = False
    for rule in rules:
        try:
            rel_to_rule = path.relative_to(rule.base).as_posix()
        except ValueError:
            continue
        if _match_rule(rule, rel_to_rule, path.is_dir()):
            ignored = not rule.negated
    return ignored


def _is_artifact_path(path: pathlib.Path) -> bool:
    rel = path.relative_to(root)
    return len(rel.parts) >= 2 and rel.parts[0] == "python" and rel.parts[1] == "artifacts"


def _is_allowed(path: pathlib.Path) -> bool:
    rel = path.relative_to(root)
    rel_posix = rel.as_posix()
    if not rel.parts:
        return False
    top = rel.parts[0]
    if top == "docs":
        if path.suffix == ".md":
            return True
        return False
    if top == "paper":
        if rel_posix.startswith("paper/build/"):
            return False
        if rel_posix.startswith("paper/preprints/Definitions/") and path.suffix in {
            ".cls",
            ".sty",
            ".bst",
            ".tex",
            ".eps",
            ".pdf",
        }:
            return True
        if path.suffix in {".tex", ".bib", ".sty", ".cls", ".sh", ".json"}:
            return True
        if path.name == "latexmkrc":
            return True
        return False
    if top == "python":
        if _is_artifact_path(path):
            # Include all produced artifacts for reproducibility, regardless of extension.
            return True
        return path.suffix in allowed_python_ext
    if top == "formal":
        return path.suffix == ".lean" or path.name in allowed_formal_names
    if top == "scripts":
        return path.suffix in allowed_scripts_ext
    return False


rules = _collect_rules(root)
files: list[pathlib.Path] = []
for base in allowed_roots:
    if not base.exists():
        continue
    for dirpath, dirnames, filenames in os.walk(base):
        dirnames[:] = [d for d in dirnames if d not in excluded_dirs]
        for filename in filenames:
            path = pathlib.Path(dirpath) / filename
            if not path.exists():
                continue
            if not _is_allowed(path):
                continue
            if not _is_artifact_path(path) and _is_ignored(path, rules):
                continue
            files.append(path)

for name in allowed_root_files:
    path = root / name
    if path.exists() and not _is_ignored(path, rules):
        files.append(path)

if not files:
    raise SystemExit("No files to package (all files ignored).")

with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
    for path in files:
        rel = path.relative_to(root).as_posix()
        zf.write(path, rel)

version_file.write_text(str(next_version), encoding="utf-8")
previous_path = None
if current_version >= 0:
    previous_path = root / f"repo_snapshot_v{current_version}.zip"
if previous_path and previous_path.exists():
    previous_path.unlink()

print(f"Wrote {zip_path}")
PY
