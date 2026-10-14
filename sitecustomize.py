import os
import sys

# Ensure src-layout packages are importable when running from repo root.
_repo_root = os.path.dirname(__file__)
_src_path = os.path.join(_repo_root, "python", "src")
if os.path.isdir(_src_path) and _src_path not in sys.path:
    sys.path.insert(0, _src_path)
