"""Fourier-space operators for periodic Navier-Stokes on T^d."""

from __future__ import annotations

import numpy as np


def kgrid_2d(N: int, L: float = 2 * np.pi) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return kx, ky, k2 on an N x N grid for a domain of length L."""
    kk = 2 * np.pi * np.fft.fftfreq(N, d=L / N)
    kx, ky = np.meshgrid(kk, kk, indexing="ij")
    k2 = kx * kx + ky * ky
    return kx, ky, k2


def kgrid_3d(N: int, L: float = 2 * np.pi) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return kx, ky, kz, k2 on an N x N x N grid for a domain of length L."""
    kk = 2 * np.pi * np.fft.fftfreq(N, d=L / N)
    kx, ky, kz = np.meshgrid(kk, kk, kk, indexing="ij")
    k2 = kx * kx + ky * ky + kz * kz
    return kx, ky, kz, k2


def dealias_mask_3d(N: int, frac: float = 2 / 3) -> np.ndarray:
    """Return a 3D de-aliasing mask using the 2/3 rule."""
    kk = np.fft.fftfreq(N)
    kcut = frac * (N // 2)
    kx, ky, kz = np.meshgrid(kk * N, kk * N, kk * N, indexing="ij")
    mask = (np.abs(kx) <= kcut) & (np.abs(ky) <= kcut) & (np.abs(kz) <= kcut)
    return mask.astype(float)


def project_div_free_2d(
    u_hat: np.ndarray,
    kx: np.ndarray,
    ky: np.ndarray,
    *,
    eps: float = 0.0,
) -> np.ndarray:
    """Project a 2D Fourier velocity field onto divergence-free modes."""
    if u_hat.shape[0] != 2:
        raise ValueError("u_hat must have shape (2, N, N)")

    k2 = kx * kx + ky * ky
    mask = k2 > eps

    u0 = u_hat[0]
    u1 = u_hat[1]

    factor = np.zeros_like(k2, dtype=complex)
    factor[mask] = (kx[mask] * u0[mask] + ky[mask] * u1[mask]) / k2[mask]

    u_proj = np.empty_like(u_hat)
    u_proj[0] = u0 - kx * factor
    u_proj[1] = u1 - ky * factor
    return u_proj


def project_div_free_3d(
    u_hat: np.ndarray,
    kx: np.ndarray,
    ky: np.ndarray,
    kz: np.ndarray,
    *,
    eps: float = 0.0,
) -> np.ndarray:
    """Project a 3D Fourier velocity field onto divergence-free modes."""
    if u_hat.shape[0] != 3:
        raise ValueError("u_hat must have shape (3, N, N, N)")

    k2 = kx * kx + ky * ky + kz * kz
    mask = k2 > eps

    u0 = u_hat[0]
    u1 = u_hat[1]
    u2 = u_hat[2]

    factor = np.zeros_like(k2, dtype=complex)
    factor[mask] = (
        kx[mask] * u0[mask] + ky[mask] * u1[mask] + kz[mask] * u2[mask]
    ) / k2[mask]

    u_proj = np.empty_like(u_hat)
    u_proj[0] = u0 - kx * factor
    u_proj[1] = u1 - ky * factor
    u_proj[2] = u2 - kz * factor
    return u_proj


def laplacian_multiplier(k2: np.ndarray) -> np.ndarray:
    """Multiplier for the Laplacian in Fourier space: Delta -> -|k|^2."""
    return -k2


def stokes_multiplier(k2: np.ndarray) -> np.ndarray:
    """Multiplier for the Stokes operator in Fourier space: A u_hat = |k|^2 P_k u_hat."""
    return k2


def hypervisc_multiplier(k2: np.ndarray, alpha: float) -> np.ndarray:
    """Multiplier for (-Delta)^alpha in Fourier space: |k|^(2 alpha)."""
    return k2**alpha
