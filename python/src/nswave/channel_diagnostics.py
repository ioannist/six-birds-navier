"""Low-rank/channel diagnostics for shell snapshot matrices."""

from __future__ import annotations

import numpy as np

Array = np.ndarray


def effective_rank_m_energy(svals: Array, energy_threshold: float) -> int:
    """
    Return smallest m whose cumulative spectral energy reaches threshold.

    Energy is measured by singular-value squares:
      sum_{i<=m} s_i^2 / sum_i s_i^2 >= energy_threshold.
    """
    s = np.asarray(svals, dtype=float).reshape(-1)
    if s.size == 0:
        return 0
    if not (0.0 < energy_threshold <= 1.0):
        raise ValueError("energy_threshold must lie in (0, 1]")
    s2 = s * s
    total = float(np.sum(s2))
    if total <= 0.0:
        return 0
    cum = np.cumsum(s2) / total
    return int(np.searchsorted(cum, energy_threshold) + 1)


def effective_rank_m95_pooled(X: Array, energy_threshold: float = 0.95) -> int:
    """
    Canonical pooled m95 from a snapshot matrix X.

    Parameters
    ----------
    X:
        Snapshot matrix with shape (n_samples, n_features).
    energy_threshold:
        Cumulative singular-value energy target.
    """
    arr = np.asarray(X, dtype=float)
    if arr.ndim != 2:
        raise ValueError("X must be a 2D matrix")
    if arr.shape[0] == 0 or arr.shape[1] == 0:
        return 0
    svals = np.linalg.svd(arr, full_matrices=False, compute_uv=False)
    return effective_rank_m_energy(svals, energy_threshold=energy_threshold)


def _fix_vector_sign(v: Array) -> Array:
    """Deterministically orient vector by largest-magnitude component sign."""
    if v.size == 0:
        return v
    idx = int(np.argmax(np.abs(v)))
    if v[idx] < 0:
        return -v
    return v


def pca_channels(X: Array, k: int) -> tuple[Array, Array]:
    """
    Compute singular values and top-k right-singular vectors.

    Parameters
    ----------
    X:
        Data matrix with shape (n_snapshots, n_features).
    k:
        Number of right-singular vectors to return.

    Returns
    -------
    svals:
        Full singular value vector.
    V:
        Top-k right-singular vectors, shape (k_eff, n_features),
        with deterministic sign orientation.
    """
    arr = np.asarray(X, dtype=float)
    if arr.ndim != 2:
        raise ValueError("X must be a 2D matrix")
    if k < 0:
        raise ValueError("k must be nonnegative")

    if arr.shape[0] == 0 or arr.shape[1] == 0 or k == 0:
        return np.zeros(0, dtype=float), np.zeros((0, arr.shape[1]), dtype=float)

    _, svals, vh = np.linalg.svd(arr, full_matrices=False)
    k_eff = min(int(k), int(vh.shape[0]))
    V = np.array(vh[:k_eff], dtype=float, copy=True)
    for i in range(k_eff):
        V[i] = _fix_vector_sign(V[i])
    return np.asarray(svals, dtype=float), V
