"""Kramers-Kronig helpers for causal, passive boundary response models."""

from __future__ import annotations

import numpy as np


def trapz_weights(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    if x.ndim != 1:
        raise ValueError("x must be a 1D array")
    if x.size < 2:
        raise ValueError("x must contain at least two points")
    dx = np.diff(x)
    if not np.all(dx > 0):
        raise ValueError("x must be strictly increasing")
    w = np.empty_like(x)
    w[0] = 0.5 * dx[0]
    w[-1] = 0.5 * dx[-1]
    if x.size > 2:
        w[1:-1] = 0.5 * (x[2:] - x[:-2])
    return w


def kk_reconstruct_re_from_im(
    omega: np.ndarray,
    imZ: np.ndarray,
    *,
    re_inf: float,
) -> np.ndarray:
    omega = np.asarray(omega, dtype=float)
    imZ = np.asarray(imZ, dtype=float)
    if omega.ndim != 1 or imZ.ndim != 1:
        raise ValueError("omega and imZ must be 1D arrays")
    if omega.size != imZ.size:
        raise ValueError("omega and imZ must have the same length")

    weights = trapz_weights(omega)
    omega_sq = omega**2
    denom = omega_sq[None, :] - omega_sq[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        integrand = (omega[None, :] * imZ[None, :]) / denom
    np.fill_diagonal(integrand, 0.0)

    re_pred = re_inf + (2.0 / np.pi) * np.sum(weights[None, :] * integrand, axis=1)
    return re_pred


def kk_static_re0_from_im(omega: np.ndarray, imZ: np.ndarray, *, re_inf: float) -> float:
    omega = np.asarray(omega, dtype=float)
    imZ = np.asarray(imZ, dtype=float)
    if omega.ndim != 1 or imZ.ndim != 1:
        raise ValueError("omega and imZ must be 1D arrays")
    if omega.size != imZ.size:
        raise ValueError("omega and imZ must have the same length")
    weights = trapz_weights(omega)
    integrand = np.zeros_like(omega)
    mask = omega != 0
    integrand[mask] = imZ[mask] / omega[mask]
    value = re_inf + (2.0 / np.pi) * np.sum(weights * integrand)
    return float(value)
