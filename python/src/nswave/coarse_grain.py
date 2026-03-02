"""Coarse-graining utilities for Fourier-space fields on T^2."""

from __future__ import annotations

import numpy as np


Array = np.ndarray


def restrict_hat_2d(omega_hat_hi: Array, N_co: int) -> Array:
    """Restrict Fourier coefficients to a smaller grid with normalization scaling."""
    N_hi = omega_hat_hi.shape[0]
    if omega_hat_hi.shape[0] != omega_hat_hi.shape[1]:
        raise ValueError("omega_hat_hi must be square")
    if N_hi % 2 != 0 or N_co % 2 != 0:
        raise ValueError("N_hi and N_co must be even")
    if N_co > N_hi:
        raise ValueError("N_co must be <= N_hi")

    hi_shift = np.fft.fftshift(omega_hat_hi)
    c0 = N_hi // 2
    h = N_co // 2
    crop = hi_shift[c0 - h : c0 + h, c0 - h : c0 + h]
    omega_hat_co = np.fft.ifftshift(crop)

    scale = (N_co / N_hi) ** 2
    return omega_hat_co * scale


def restrict_hat_3d(arr_hat_hi: Array, N_co: int) -> Array:
    """Restrict 3D Fourier coefficients with normalization scaling."""
    if arr_hat_hi.ndim == 3:
        spatial_axes = (0, 1, 2)
        N_hi = arr_hat_hi.shape[0]
        if arr_hat_hi.shape[1] != N_hi or arr_hat_hi.shape[2] != N_hi:
            raise ValueError("arr_hat_hi must be cubic for scalar input")
    elif arr_hat_hi.ndim == 4:
        spatial_axes = (1, 2, 3)
        N_hi = arr_hat_hi.shape[1]
        if arr_hat_hi.shape[2] != N_hi or arr_hat_hi.shape[3] != N_hi:
            raise ValueError("arr_hat_hi must be cubic for vector input")
    else:
        raise ValueError("arr_hat_hi must have ndim 3 or 4")

    if N_hi % 2 != 0 or N_co % 2 != 0:
        raise ValueError("N_hi and N_co must be even")
    if N_co > N_hi:
        raise ValueError("N_co must be <= N_hi")

    hi_shift = np.fft.fftshift(arr_hat_hi, axes=spatial_axes)
    c0 = N_hi // 2
    h = N_co // 2
    if arr_hat_hi.ndim == 3:
        crop = hi_shift[c0 - h : c0 + h, c0 - h : c0 + h, c0 - h : c0 + h]
    else:
        crop = hi_shift[:, c0 - h : c0 + h, c0 - h : c0 + h, c0 - h : c0 + h]
    arr_hat_co = np.fft.ifftshift(crop, axes=spatial_axes)

    scale = (N_co / N_hi) ** 3
    return arr_hat_co * scale


def shell_index_from_k2(k2: Array) -> Array:
    """Return integer shell index s = floor(sqrt(k2)+0.5)."""
    return np.floor(np.sqrt(k2) + 0.5).astype(int)
