import os
import sys

_root = os.path.dirname(os.path.dirname(__file__))
_src_path = os.path.join(_root, "src")
if os.path.isdir(_src_path) and _src_path not in sys.path:
    sys.path.insert(0, _src_path)
