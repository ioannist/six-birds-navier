import os
import sys
from pkgutil import extend_path

# Shim to make src-layout package importable from repo root.
_repo_root = os.path.dirname(os.path.dirname(__file__))
_src_path = os.path.join(_repo_root, "python", "src")
if _src_path not in sys.path:
    sys.path.insert(0, _src_path)

__path__ = extend_path(__path__, __name__)
