"""Localization diagnostics for real-space 3D vector fields."""

from __future__ import annotations

import numpy as np

Array = np.ndarray


def max_local_energy_fraction(u: Array, window: int) -> float:
    """
    Return the maximum periodic cube-local energy fraction.

    Parameters
    ----------
    u:
        Real-space vector field with shape (3, N, N, N).
    window:
        Cube side length in grid cells.
    """
    if u.ndim != 4 or u.shape[0] != 3:
        raise ValueError("u must have shape (3, N, N, N)")
    if u.shape[1] != u.shape[2] or u.shape[1] != u.shape[3]:
        raise ValueError("u must be cubic in space")
    if window < 1:
        raise ValueError("window must be >= 1")

    n = u.shape[1]
    w = min(int(window), n)

    e = 0.5 * np.sum(u * u, axis=0)
    e_total = float(np.sum(e))
    if e_total <= 0.0:
        return 0.0

    kernel = np.zeros((n, n, n), dtype=float)
    kernel[:w, :w, :w] = 1.0
    local_sum = np.fft.ifftn(np.fft.fftn(e) * np.fft.fftn(kernel)).real
    max_local = float(np.max(local_sum))
    frac = max_local / e_total
    return float(min(max(frac, 0.0), 1.0))


def baseline_volume_fraction(N: int, window: int) -> float:
    """Return cube-volume fraction (window^3 / N^3) with clipping to [1, N]."""
    if N <= 0:
        raise ValueError("N must be positive")
    if window < 1:
        raise ValueError("window must be >= 1")
    w = min(int(window), N)
    return float((w**3) / (N**3))


def localization_factors(
    conc: Array,
    *,
    N: int,
    windows: Array,
    j_values: Array,
) -> tuple[Array, Array, Array]:
    """
    Compute baseline fraction and normalized localization factors.

    Parameters
    ----------
    conc:
        Concentration values with shell index on the last axis.
    N:
        Grid size.
    windows:
        Window sizes per shell.
    j_values:
        Shell indices per shell.

    Returns
    -------
    f:
        Baseline volume fractions per shell.
    L:
        Excess localization factor conc / f.
    kappa:
        Channelized scaling factor L / (j+1).
    """
    c = np.asarray(conc, dtype=float)
    w = np.asarray(windows, dtype=float)
    j = np.asarray(j_values, dtype=float)
    if c.shape[-1] != w.shape[0] or w.shape != j.shape:
        raise ValueError("conc last axis must match windows/j_values length")
    if N <= 0:
        raise ValueError("N must be positive")

    f = (w**3) / float(N**3)
    L = c / np.maximum(f, 1e-30)
    kappa = L / (j + 1.0)
    return f, L, kappa


def window_for_shell(N: int, j: int, c: float = 1.0) -> int:
    """
    Return cube window side length for shell scale ~2^{-j}.

    Uses round(c * N / 2**j), clipped to [1, N].
    """
    if N <= 0:
        raise ValueError("N must be positive")
    if c <= 0:
        raise ValueError("c must be positive")
    w = int(round(c * N / (2.0 ** j)))
    return max(1, min(N, w))
