#!/usr/bin/env bash
set -euo pipefail

MODE="--fast"
INCLUDE_SLOW="0"

for arg in "$@"; do
  case "$arg" in
    --include-slow)
      INCLUDE_SLOW="1"
      ;;
    --fast)
      MODE="--fast"
      ;;
    --no-figures)
      MODE="--no-figures"
      ;;
    *)
      echo "Unknown arg: $arg"
      echo "Usage: build_paper.sh [--include-slow] [--fast] [--no-figures]"
      exit 2
      ;;
  esac
done

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

# Figures (paper artifacts)
if [ "$MODE" != "--no-figures" ]; then
  if [ "$INCLUDE_SLOW" = "1" ]; then
    python python/scripts/ns_paper_figures.py --include-slow --clean
  else
    python python/scripts/ns_paper_figures.py --clean
  fi
fi

# Build LaTeX
cd paper
latexmk -pdf -quiet main.tex

echo "Built PDF: paper/build/main.pdf"
