"""Dyadic shell utilities for 3D Fourier fields (diagnostic layer)."""

from __future__ import annotations

import numpy as np

Array = np.ndarray


def dyadic_shell_index(k_mag: Array, *, k0: float = 1.0) -> Array:
    """Return dyadic shell index for each k (k=0 mapped to -1)."""
    idx = -np.ones_like(k_mag, dtype=int)
    mask = k_mag >= k0
    idx[mask] = np.floor(np.log2(k_mag[mask] / k0)).astype(int)
    return idx


def dyadic_shell_mask(k_mag: Array, j: int, *, k0: float = 1.0) -> Array:
    """Return mask for shell j: k in [k0*2^j, k0*2^(j+1))."""
    k_lo = k0 * (2.0 ** j)
    k_hi = k0 * (2.0 ** (j + 1))
    return (k_mag >= k_lo) & (k_mag < k_hi)


def _energy_scale(u_hat: Array, *, L: float) -> float:
    """Return Parseval scaling for energy with numpy FFT conventions."""
    if u_hat.ndim < 3:
        raise ValueError("u_hat must have at least 3 spatial dims")
    if u_hat.ndim in (3, 4) and u_hat.shape[0] in (2, 3):
        d = u_hat.ndim - 1
    else:
        d = u_hat.ndim
    N = u_hat.shape[-1]
    return L**d / (N ** (2 * d))


def shell_energy(u_hat: Array, shell_mask: Array, *, L: float) -> float:
    """Energy in a shell using Parseval scaling (diagnostic; real-space convention)."""
    if u_hat.ndim == shell_mask.ndim + 1:
        masked = u_hat * shell_mask[None, ...]
    else:
        masked = u_hat * shell_mask
    scale = _energy_scale(u_hat, L=L)
    return 0.5 * float(np.sum(np.abs(masked) ** 2) * scale)


def shell_energies(u_hat: Array, k_mag: Array, *, j_max: int | None = None, L: float = 2 * np.pi, k0: float = 1.0) -> Array:
    """Return dyadic shell energies E_j for j=0..j_max (k=0 excluded)."""
    k_max = float(np.max(k_mag))
    if j_max is None:
        if k_max < k0:
            j_max = 0
        else:
            j_max = int(np.floor(np.log2(k_max / k0)))
    energies = np.zeros(j_max + 1)
    for j in range(j_max + 1):
        mask = dyadic_shell_mask(k_mag, j, k0=k0)
        energies[j] = shell_energy(u_hat, mask, L=L)
    return energies


def frontier_scale(E_shell: Array, threshold: float) -> int:
    """Return max j with E_j >= threshold, else -1."""
    idx = np.where(E_shell >= threshold)[0]
    return int(idx.max()) if idx.size else -1


def percentile_frontier(E_shell: Array, p: float) -> int:
    """Return smallest j with tail-energy fraction <= p (p in [0,1])."""
    total = float(np.sum(E_shell))
    if total <= 0:
        return -1
    tail_from_j = np.cumsum(E_shell[::-1])[::-1]
    frac = tail_from_j / total
    idx = np.where(frac <= p)[0]
    if idx.size == 0:
        return int(len(E_shell) - 1)
    return int(idx[0])


def diagnostic_transfer_proxy(u_hat: Array, nl_hat: Array, k_mag: Array, j: int, *, L: float = 2 * np.pi, k0: float = 1.0) -> float:
    """Diagnostic cumulative transfer proxy across shells <= j (not a proof object)."""
    k_hi = k0 * (2.0 ** (j + 1))
    mask = (k_mag >= k0) & (k_mag < k_hi)
    if u_hat.ndim == mask.ndim + 1:
        dot = np.sum(np.conj(u_hat) * nl_hat, axis=0)
    else:
        dot = np.conj(u_hat) * nl_hat
    scale = _energy_scale(u_hat, L=L)
    return float(np.real(np.sum(dot * mask) * scale))


def monotone_envelope(j_series: Array) -> Array:
    """Return the monotone envelope of a frontier sequence (nondecreasing)."""
    return np.maximum.accumulate(j_series)


def crossing_times(t: Array, j_mon: Array, *, j_start: int | None = None, j_end: int | None = None) -> tuple[Array, Array]:
    """Return first-hit times for each integer j in [j_start, j_end]."""
    if t.shape != j_mon.shape:
        raise ValueError("t and j_mon must have the same shape")
    j0 = int(j_mon[0]) if j_start is None else int(j_start)
    jmax = int(np.max(j_mon)) if j_end is None else int(j_end)
    j_vals = np.arange(j0, jmax + 1, dtype=int)
    t_reach = np.full_like(j_vals, np.nan, dtype=float)
    for idx, j in enumerate(j_vals):
        hit = np.where(j_mon >= j)[0]
        if hit.size:
            t_reach[idx] = float(t[hit[0]])
    return j_vals, t_reach
