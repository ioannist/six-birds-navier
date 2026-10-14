#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PAPER_DIR="$ROOT_DIR/paper"
BUILD_DIR="$PAPER_DIR/build"
STAGE_DIR="$BUILD_DIR/preprints_source_staging"
ZIP_PATH="$BUILD_DIR/preprints_source_upload.zip"
PREVIEW_PDF="$BUILD_DIR/preprints_preview.pdf"
PREPRINTS_MAIN_PDF="$BUILD_DIR/main_preprints.pdf"
ARCHIVE_ROOT="$PAPER_DIR/archive"
TEMPLATE_DIR="$PAPER_DIR/preprints"
DEFINITIONS_DIR="$TEMPLATE_DIR/Definitions"
PREAMBLE_FILE="$TEMPLATE_DIR/preamble.tex"
MAIN_TEX="$PAPER_DIR/main.tex"
REFS_BIB="$PAPER_DIR/refs.bib"
SECTIONS_DIR="$PAPER_DIR/sections"

if [[ ! -f "$MAIN_TEX" ]]; then
  echo "[package_preprints] ERROR: missing $MAIN_TEX" >&2
  exit 2
fi
if [[ ! -f "$REFS_BIB" ]]; then
  echo "[package_preprints] ERROR: missing $REFS_BIB" >&2
  exit 2
fi
if [[ ! -d "$SECTIONS_DIR" ]]; then
  echo "[package_preprints] ERROR: missing $SECTIONS_DIR" >&2
  exit 2
fi
if [[ ! -d "$DEFINITIONS_DIR" ]]; then
  echo "[package_preprints] ERROR: missing template definitions at $DEFINITIONS_DIR" >&2
  exit 2
fi
if [[ ! -f "$PREAMBLE_FILE" ]]; then
  echo "[package_preprints] ERROR: missing template preamble at $PREAMBLE_FILE" >&2
  exit 2
fi

# Keep appendix figure table in sync with manifest before packaging.
python "$ROOT_DIR/python/scripts/ns_paper_appendix_figure_table.py"

mkdir -p "$ARCHIVE_ROOT"
ARCHIVE_TAG="${1:-}"
if [[ -z "$ARCHIVE_TAG" ]]; then
  max=0
  shopt -s nullglob
  for dir in "$ARCHIVE_ROOT"/preprints_v*; do
    base="$(basename "$dir")"
    if [[ "$base" =~ ^preprints_v([0-9]+)$ ]]; then
      n="${BASH_REMATCH[1]}"
      if (( n > max )); then
        max="$n"
      fi
    fi
  done
  shopt -u nullglob
  ARCHIVE_TAG="preprints_v$((max + 1))"
fi

ARCHIVE_DIR="$ARCHIVE_ROOT/$ARCHIVE_TAG"
rm -rf "$ARCHIVE_DIR"
mkdir -p "$ARCHIVE_DIR/current_source/sections" "$ARCHIVE_DIR/current_source/figures" "$ARCHIVE_DIR/preprints_bundle"

rm -rf "$BUILD_DIR"
mkdir -p "$STAGE_DIR/sections" "$STAGE_DIR/figures"

cp "$SECTIONS_DIR/"*.tex "$STAGE_DIR/sections/"
cp "$REFS_BIB" "$STAGE_DIR/refs.bib"
cp "$PREAMBLE_FILE" "$STAGE_DIR/preamble.tex"
cp -r "$DEFINITIONS_DIR" "$STAGE_DIR/Definitions"

python - "$ROOT_DIR" "$PAPER_DIR" "$STAGE_DIR" <<'PY'
from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

root = Path(sys.argv[1])
paper_dir = Path(sys.argv[2])
stage_dir = Path(sys.argv[3])
stage_sections = stage_dir / "sections"
stage_figures = stage_dir / "figures"

main_tex = (paper_dir / "main.tex").read_text(encoding="utf-8", errors="replace")
abstract_tex = (paper_dir / "sections" / "00_title_abstract.tex").read_text(
    encoding="utf-8", errors="replace"
)

title_m = re.search(r"\\title\{(.*?)\}", main_tex, re.S)
if not title_m:
    raise SystemExit("[package_preprints] ERROR: could not extract title from paper/main.tex")
title = " ".join(title_m.group(1).split())

abstract_m = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", abstract_tex, re.S)
if not abstract_m:
    raise SystemExit(
        "[package_preprints] ERROR: could not extract abstract from paper/sections/00_title_abstract.tex"
    )
abstract = " ".join(abstract_m.group(1).split())

inputs = re.findall(r"\\input\{(sections/[^}]+)\}", main_tex)
if not inputs:
    raise SystemExit("[package_preprints] ERROR: no section inputs found in paper/main.tex")
inputs = [s for s in inputs if s != "sections/00_title_abstract.tex" and s != "sections/00_title_abstract"]

include_re = re.compile(r"(\\includegraphics(?:\[[^\]]*\])?\{)([^}]+)(\})")
source_for_basename: dict[str, Path] = {}
missing: list[Path] = []

for section_path in stage_sections.glob("*.tex"):
    text = section_path.read_text(encoding="utf-8", errors="replace")

    def repl(match: re.Match[str]) -> str:
        prefix, path_str, suffix = match.groups()
        src = (paper_dir / path_str).resolve()
        if not src.exists() or not src.is_file():
            missing.append(paper_dir / path_str)
            return match.group(0)

        name = src.name
        prev = source_for_basename.get(name)
        if prev is not None and prev != src:
            raise SystemExit(
                f"[package_preprints] ERROR: figure basename collision: {prev} and {src}"
            )
        source_for_basename[name] = src
        dst = stage_figures / name
        if not dst.exists():
            shutil.copy2(src, dst)
        return f"{prefix}figures/{name}{suffix}"

    text = include_re.sub(repl, text)
    section_path.write_text(text, encoding="utf-8")

if missing:
    first = missing[0]
    raise SystemExit(f"[package_preprints] ERROR: missing referenced figure file: {first}")

keywords = (
    "navier--stokes; no-zeno; six birds; anti-localization; route capacity; "
    "mechanized theorem proving; lean"
)

content = [
    r"\documentclass[apajournal,article,submit,moreauthors,pdftex]{Definitions/mdpi}",
    r"\input{preamble}",
    "",
    r"\firstpage{1}",
    r"\makeatletter",
    r"\setcounter{page}{\@firstpage}",
    r"\makeatother",
    r"\pubvolume{1}",
    r"\issuenum{1}",
    r"\articlenumber{0}",
    r"\pubyear{2026}",
    r"\copyrightyear{2026}",
    r"\datereceived{ }",
    r"\daterevised{ }",
    r"\dateaccepted{ }",
    r"\datepublished{ }",
    "",
    f"\\Title{{{title}}}",
    r"\newcommand{\orcidauthorA}{0009-0009-7659-5964}",
    r"\Author{Ioannis Tsiokos $^{1}$\orcidA{}}",
    r"\AuthorNames{Ioannis Tsiokos}",
    r"\address{$^{1}$ \quad Automorph Inc., Wilmington, DE 19806, USA}",
    r"\corres{Correspondence: ioannis@automorph.io}",
    f"\\abstract{{{abstract}}}",
    f"\\keyword{{{keywords}}}",
    "",
    r"\begin{document}",
    "",
]

for inp in inputs:
    content.append(f"\\input{{{inp.removesuffix('.tex')}}}")

content += [
    "",
    r"\authorcontributions{I.T.\ conceived the study, implemented experiments and formal verification, and wrote the manuscript.}",
    r"\funding{This research received no external funding.}",
    r"\institutionalreview{Not applicable.}",
    r"\informedconsent{Not applicable.}",
    r"\dataavailability{All data, code, and reproducibility artifacts are available in the companion repository.}",
    r"\useofartificialintelligence{The author used AI-assisted tools for coding and editorial support during repository development and manuscript preparation. All technical claims, theorem statements, proofs, numerical analyses, and final wording were reviewed and validated by the author, who takes full responsibility for the content.}",
    r"\acknowledgments{The author acknowledges the use of Lean 4/mathlib and standard Python scientific libraries.}",
    r"\conflictsofinterest{The author declares no conflicts of interest.}",
    "",
    r"\begin{adjustwidth}{-\extralength}{0cm}",
    r"\reftitle{References}",
    r"\externalbibliography{yes}",
    r"\bibliography{refs}",
    r"\PublishersNote{}",
    r"\end{adjustwidth}",
    "",
    r"\end{document}",
    "",
]

(stage_dir / "main.tex").write_text("\n".join(content), encoding="utf-8")
PY

if command -v latexmk >/dev/null 2>&1; then
  (cd "$STAGE_DIR" && latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex >/dev/null)
elif command -v pdflatex >/dev/null 2>&1; then
  (cd "$STAGE_DIR" && pdflatex -interaction=nonstopmode -halt-on-error main.tex >/dev/null)
  (cd "$STAGE_DIR" && bibtex main >/dev/null)
  (cd "$STAGE_DIR" && pdflatex -interaction=nonstopmode -halt-on-error main.tex >/dev/null)
  (cd "$STAGE_DIR" && pdflatex -interaction=nonstopmode -halt-on-error main.tex >/dev/null)
else
  echo "[package_preprints] ERROR: no LaTeX builder found (latexmk/pdflatex)." >&2
  exit 3
fi

if [[ ! -f "$STAGE_DIR/main.bbl" ]]; then
  echo "[package_preprints] ERROR: template build did not produce main.bbl" >&2
  exit 4
fi
if [[ ! -f "$STAGE_DIR/main.pdf" ]]; then
  echo "[package_preprints] ERROR: template build did not produce main.pdf" >&2
  exit 4
fi

cp "$STAGE_DIR/main.pdf" "$PREVIEW_PDF"
cp "$STAGE_DIR/main.pdf" "$PREPRINTS_MAIN_PDF"

# Archive current canonical paper source snapshot.
cp "$MAIN_TEX" "$ARCHIVE_DIR/current_source/main.tex"
cp "$PAPER_DIR/macros.tex" "$ARCHIVE_DIR/current_source/macros.tex"
cp "$REFS_BIB" "$ARCHIVE_DIR/current_source/refs.bib"
cp "$PAPER_DIR/figures_manifest.json" "$ARCHIVE_DIR/current_source/figures_manifest.json"
cp "$SECTIONS_DIR/"*.tex "$ARCHIVE_DIR/current_source/sections/"
cp -r "$TEMPLATE_DIR" "$ARCHIVE_DIR/current_source/preprints"
cp "$STAGE_DIR/figures/"* "$ARCHIVE_DIR/current_source/figures/" 2>/dev/null || true

find "$STAGE_DIR" -name "*.aux" -delete
find "$STAGE_DIR" -name "*.log" -delete
find "$STAGE_DIR" -name "*.out" -delete
find "$STAGE_DIR" -name "*.toc" -delete
find "$STAGE_DIR" -name "*.fdb_latexmk" -delete
find "$STAGE_DIR" -name "*.fls" -delete
rm -f "$STAGE_DIR/main.pdf" "$STAGE_DIR/main.blg"
find "$STAGE_DIR" -name ".DS_Store" -delete

(
  cd "$STAGE_DIR"
  zip -r "$ZIP_PATH" . >/dev/null
)

cp -r "$STAGE_DIR" "$ARCHIVE_DIR/preprints_bundle/source"
cp "$ZIP_PATH" "$ARCHIVE_DIR/preprints_bundle/preprints_source_upload.zip"
cp "$PREVIEW_PDF" "$ARCHIVE_DIR/preprints_bundle/preprints_preview.pdf"
cp "$PREPRINTS_MAIN_PDF" "$ARCHIVE_DIR/preprints_bundle/main_preprints.pdf"

echo "[package_preprints] Wrote $ZIP_PATH"
echo "[package_preprints] Wrote preview PDF $PREVIEW_PDF"
echo "[package_preprints] Wrote canonical preprints PDF $PREPRINTS_MAIN_PDF"
echo "[package_preprints] Archived current source at $ARCHIVE_DIR/current_source"
echo "[package_preprints] Archived preprints bundle at $ARCHIVE_DIR/preprints_bundle"
