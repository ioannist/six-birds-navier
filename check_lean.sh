#!/usr/bin/env bash
set -euo pipefail

(cd formal && lake build)
