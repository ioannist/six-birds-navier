"""Inference helpers for passive Debye boundary models."""

from __future__ import annotations

from typing import Tuple

import numpy as np
from scipy.optimize import least_squares


def Z_debye_k1(omega: np.ndarray, a0: float, a1: float, b1: float) -> np.ndarray:
    omega = np.asarray(omega, dtype=float)
    return a0 + a1 / (1.0 - 1j * omega / b1)


def R_pred_grid(
    omega: np.ndarray,
    R0: np.ndarray,
    r: np.ndarray,
    t: np.ndarray,
    Dx: float,
) -> np.ndarray:
    omega = np.asarray(omega, dtype=float)
    phase = np.exp(2j * omega * Dx)
    return r * np.exp(-2j * omega * Dx) + (t * t * R0) / (1.0 - R0 * r * phase)


def gaussian_source(omega: np.ndarray, omega0: float, sigma: float) -> np.ndarray:
    omega = np.asarray(omega, dtype=float)
    return np.exp(-0.5 * ((omega - omega0) / sigma) ** 2)


def synthesize_waveform_irfft(
    *,
    N: int,
    dt: float,
    S_pos: np.ndarray,
    Rpred_pos: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    omega_pos = 2.0 * np.pi * np.fft.rfftfreq(N, dt)
    F_pos = np.asarray(S_pos, dtype=complex) * np.asarray(Rpred_pos, dtype=complex)
    F_pos = F_pos.copy()
    if F_pos.size > 0:
        F_pos[0] = 0.0
    y = np.fft.irfft(F_pos, n=N)
    t = np.arange(N, dtype=float) * dt
    return t, y, F_pos


def estimate_Rpred_from_waveform(
    *,
    waveform: np.ndarray,
    S_pos: np.ndarray,
    eps: float = 1e-12,
) -> np.ndarray:
    F_meas = np.fft.rfft(np.asarray(waveform, dtype=float))
    S_pos = np.asarray(S_pos, dtype=complex)
    mask = np.abs(S_pos) > 1e-3 * np.max(np.abs(S_pos))
    R_est = np.zeros_like(F_meas)
    R_est[mask] = F_meas[mask] / (S_pos[mask] + eps)
    return R_est


def estimate_Rpred_from_baseline_ratio(
    *,
    waveform_true: np.ndarray,
    waveform_base: np.ndarray,
    R_base: np.ndarray,
    eps: float = 1e-12,
    min_base_frac: float = 1e-3,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Estimate system reflection R_true(omega) from two real waveforms with identical source:
      F_true ~= S * R_true, F_base ~= S * R_base => R_true ~= (F_true/F_base) * R_base.
    Returns (R_est, mask) where mask selects bins with sufficiently large |F_base|.
    """
    F_true = np.fft.rfft(np.asarray(waveform_true, dtype=float))
    F_base = np.fft.rfft(np.asarray(waveform_base, dtype=float))
    R_base = np.asarray(R_base, dtype=complex)
    if R_base.shape != F_true.shape:
        raise ValueError("R_base must have the same shape as rfft(waveform)")
    threshold = min_base_frac * np.max(np.abs(F_base)) if F_base.size else 0.0
    mask = np.abs(F_base) > threshold
    R_est = np.zeros_like(F_true)
    R_est[mask] = (F_true[mask] / (F_base[mask] + eps)) * R_base[mask]
    return R_est, mask


def estimate_Rout_from_two_point_spectra(
    *,
    omega_pos: np.ndarray,
    Psi1: np.ndarray,
    Psi2: np.ndarray,
    y1: float,
    y2: float,
    min_amp_frac: float = 1e-3,
    eps: float = 1e-12,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Solve for A_in, A_out from two-point plane-wave data and return R_out = A_out / A_in.
    """
    omega_pos = np.asarray(omega_pos, dtype=float)
    Psi1 = np.asarray(Psi1, dtype=complex)
    Psi2 = np.asarray(Psi2, dtype=complex)
    if omega_pos.shape != Psi1.shape or omega_pos.shape != Psi2.shape:
        raise ValueError("omega_pos, Psi1, and Psi2 must have the same shape")

    E1m = np.exp(-1j * omega_pos * y1)
    E1p = np.exp(1j * omega_pos * y1)
    E2m = np.exp(-1j * omega_pos * y2)
    E2p = np.exp(1j * omega_pos * y2)
    denom = E1m * E2p - E2m * E1p
    A_in = (Psi1 * E2p - Psi2 * E1p) / (denom + eps)
    A_out = (Psi2 * E1m - Psi1 * E2m) / (denom + eps)
    R_out = A_out / (A_in + eps)

    amp = np.abs(A_in)
    band_mask = (omega_pos >= 0.5) & (omega_pos <= 5.0)
    if np.any(band_mask):
        max_amp = np.max(amp[band_mask])
    else:
        max_amp = np.max(amp) if amp.size else 0.0
    threshold = min_amp_frac * max_amp
    mask = band_mask & (amp > threshold)
    R_out_masked = np.zeros_like(R_out)
    R_out_masked[mask] = R_out[mask]
    return R_out_masked, mask


def fit_debye_k1_from_Rpred(
    *,
    omega_pos: np.ndarray,
    Rpred_target: np.ndarray,
    r: np.ndarray,
    t: np.ndarray,
    Dx: float,
    p0: tuple[float, float, float],
    weight: np.ndarray | None = None,
    max_nfev: int = 200,
) -> dict:
    omega_pos = np.asarray(omega_pos, dtype=float)
    Rpred_target = np.asarray(Rpred_target, dtype=complex)
    r = np.asarray(r, dtype=complex)
    t = np.asarray(t, dtype=complex)

    if any(val <= 0 for val in p0):
        raise ValueError("initial guess p0 must be positive")

    if weight is None:
        w = np.ones_like(omega_pos, dtype=float)
    else:
        w = np.asarray(weight, dtype=float)
        if w.shape != omega_pos.shape:
            raise ValueError("weight must match omega_pos shape")

    def residual(u: np.ndarray) -> np.ndarray:
        a0 = float(np.exp(u[0]))
        a1 = float(np.exp(u[1]))
        b1 = float(np.exp(u[2]))
        Z = Z_debye_k1(omega_pos, a0, a1, b1)
        R0 = (1.0 - Z) / (1.0 + Z)
        R_model = R_pred_grid(omega_pos, R0, r, t, Dx)
        diff = w * (R_model - Rpred_target)
        return np.concatenate([diff.real, diff.imag])

    u0 = np.log(np.array(p0, dtype=float))
    sol = least_squares(residual, u0, max_nfev=max_nfev)
    a0 = float(np.exp(sol.x[0]))
    a1 = float(np.exp(sol.x[1]))
    b1 = float(np.exp(sol.x[2]))

    return {
        "a0": a0,
        "a1": a1,
        "b1": b1,
        "success": bool(sol.success),
        "cost": float(sol.cost),
        "nfev": int(sol.nfev),
        "message": str(sol.message),
    }


def passivity_residuals(
    omega: np.ndarray,
    Z: np.ndarray,
    R0: np.ndarray,
) -> tuple[float, float]:
    resid_R = max(0.0, float(np.max(np.abs(R0) - 1.0)))
    resid_Z = max(0.0, float(-np.min(np.real(Z))))
    return resid_R, resid_Z
