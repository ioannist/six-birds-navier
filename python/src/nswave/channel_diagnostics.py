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
    if not (0.0 < energy_threshold <= 1.0):
        raise ValueError("energy_threshold must lie in (0, 1]")
    if not np.all(np.isfinite(s)) or np.any(s < 0):
        raise ValueError("singular values must be finite and nonnegative")
    if s.size == 0 or np.max(s) == 0:
        return 0
    s = np.sort(s)[::-1] / np.max(s)
    s2 = s * s
    total = float(np.sum(s2))
    if total <= 0.0:
        return 0
    cum = np.cumsum(s2) / total
    cum[-1] = 1.0
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
    if not np.all(np.isfinite(arr)):
        raise ValueError("X must be finite")
    if arr.shape[0] == 0 or arr.shape[1] == 0:
        return 0
    scale = float(np.max(np.abs(arr)))
    if scale == 0:
        return 0
    svals = np.linalg.svd(arr / scale, full_matrices=False, compute_uv=False)
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


def channel_span_snapshot_bounds(X: Array, channels: Array, eta: float) -> tuple[Array, Array]:
    """Charge discarded energy when returning a channel-span bound to snapshots.

    Channel rows must be orthonormal. eta is a supplied upper bound on local
    energy fraction for their whole span. Returns per-snapshot residual
    fraction and full local fraction bound; both are zero on a zero snapshot.
    This does not infer exact feasibility from an approximate decomposition.
    """
    arr = np.asarray(X, dtype=float)
    basis = np.asarray(channels, dtype=float)
    if arr.ndim != 2 or basis.ndim != 2 or arr.shape[1] != basis.shape[1]:
        raise ValueError("X and channels must be matrices with matching features")
    if not np.all(np.isfinite(arr)) or not np.all(np.isfinite(basis)):
        raise ValueError("X and channels must be finite")
    if not np.isfinite(eta) or eta < 0:
        raise ValueError("eta must be finite and nonnegative")
    if not np.allclose(basis @ basis.T, np.eye(basis.shape[0]), rtol=1e-12, atol=1e-12):
        raise ValueError("channel rows must be orthonormal")
    delta = np.zeros(arr.shape[0])
    bounds = np.zeros(arr.shape[0])
    if arr.size == 0:
        return delta, bounds
    scale = float(np.max(np.abs(arr)))
    if scale == 0:
        return delta, bounds
    scaled = arr / scale
    projected = (scaled @ basis.T) @ basis
    residual = scaled - projected
    total_sq = np.sum(scaled * scaled, axis=1)
    nonzero = total_sq > 0
    delta[nonzero] = np.sum(residual[nonzero] ** 2, axis=1) / total_sq[nonzero]
    retained = np.sum(projected[nonzero] ** 2, axis=1) / total_sq[nonzero]
    bounds[nonzero] = np.minimum(
        1.0, (np.sqrt(min(eta, 1.0) * retained) + np.sqrt(delta[nonzero])) ** 2
    )
    return delta, bounds
