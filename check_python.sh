#!/usr/bin/env bash
set -euo pipefail

PYTHONPATH=python/src python3 -m pytest python/tests -q -m "not slow"
