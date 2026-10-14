"""Deterministic NPZ writer utilities."""

from __future__ import annotations

import hashlib
import io
from pathlib import Path
import zipfile

import numpy as np


def savez_deterministic(path: str | Path, **arrays: np.ndarray) -> None:
    """
    Write a deterministic .npz archive.

    Determinism controls:
    - sorted keys
    - fixed zip timestamp
    - ZIP_STORED (no compression)
    - np.save(..., allow_pickle=False)
    """
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with zipfile.ZipFile(out_path, mode="w", compression=zipfile.ZIP_STORED) as zf:
        for key in sorted(arrays.keys()):
            arr = np.asarray(arrays[key])
            bio = io.BytesIO()
            np.save(bio, arr, allow_pickle=False)
            payload = bio.getvalue()

            info = zipfile.ZipInfo(
                filename=f"{key}.npy",
                date_time=(1980, 1, 1, 0, 0, 0),
            )
            info.compress_type = zipfile.ZIP_STORED
            info.create_system = 3
            info.external_attr = 0o600 << 16
            zf.writestr(info, payload)


def sha256_file(path: str | Path) -> str:
    """Return SHA-256 hex digest for a file."""
    p = Path(path)
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()
