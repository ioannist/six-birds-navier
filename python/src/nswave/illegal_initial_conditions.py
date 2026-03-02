"""Construct intentionally SPT-illegal, wavepacket-like initial conditions."""

from __future__ import annotations

import numpy as np

from nswave.dyadic import dyadic_shell_mask
from nswave.operators import dealias_mask_3d, kgrid_3d, project_div_free_3d


def _wrap_displacement(x: np.ndarray, c: float) -> np.ndarray:
    """Periodic displacement in unit-torus coordinates mapped to [-0.5, 0.5)."""
    return ((x - c + 0.5) % 1.0) - 0.5


def _energy_from_u_hat(u_hat: np.ndarray, L: float = 2 * np.pi) -> float:
    u0 = np.fft.ifftn(u_hat[0]).real
    u1 = np.fft.ifftn(u_hat[1]).real
    u2 = np.fft.ifftn(u_hat[2]).real
    return float(0.5 * np.mean(u0 * u0 + u1 * u1 + u2 * u2) * (L**3))


def make_localized_blob_u0_hat(
    N: int,
    center_xyz: tuple[float, float, float],
    sigma: float,
    seed: int,
) -> np.ndarray:
    """
    Build a localized divergence-free initial field in Fourier space.

    Construction:
    - Define periodic Gaussian scalar blob phi(x) on unit torus.
    - Set vector potential A(x) = phi(x) * e with seeded random unit vector e.
    - Set u_hat = i k x A_hat (spectral curl), which is divergence-free.

    Parameters
    ----------
    center_xyz:
        Blob center in unit-torus coordinates [0, 1)^3.
    sigma:
        Gaussian width in unit-torus coordinates.
    """
    if N <= 0:
        raise ValueError("N must be positive")
    if sigma <= 0:
        raise ValueError("sigma must be positive")
    if len(center_xyz) != 3:
        raise ValueError("center_xyz must have length 3")

    cx, cy, cz = (float(center_xyz[0]), float(center_xyz[1]), float(center_xyz[2]))
    grid = np.arange(N, dtype=float) / float(N)
    x, y, z = np.meshgrid(grid, grid, grid, indexing="ij")
    dx = _wrap_displacement(x, cx)
    dy = _wrap_displacement(y, cy)
    dz = _wrap_displacement(z, cz)
    r2 = dx * dx + dy * dy + dz * dz
    phi = np.exp(-0.5 * r2 / (sigma * sigma))

    rng = np.random.default_rng(seed)
    e = rng.standard_normal(3)
    e_norm = float(np.linalg.norm(e))
    if e_norm <= 1e-14:
        e = np.array([1.0, 0.0, 0.0], dtype=float)
    else:
        e = e / e_norm

    A = np.empty((3, N, N, N), dtype=float)
    A[0] = e[0] * phi
    A[1] = e[1] * phi
    A[2] = e[2] * phi

    A_hat = np.empty_like(A, dtype=complex)
    for i in range(3):
        A_hat[i] = np.fft.fftn(A[i])

    kx, ky, kz, _ = kgrid_3d(N, L=2 * np.pi)
    u_hat = np.empty_like(A_hat)
    u_hat[0] = 1j * (ky * A_hat[2] - kz * A_hat[1])
    u_hat[1] = 1j * (kz * A_hat[0] - kx * A_hat[2])
    u_hat[2] = 1j * (kx * A_hat[1] - ky * A_hat[0])
    u_hat[:, 0, 0, 0] = 0.0
    return u_hat


def bandlimit_to_shell(u0_hat: np.ndarray, j: int, dealias: bool = True) -> np.ndarray:
    """
    Restrict field to dyadic shell S_j and enforce discrete divergence-free structure.
    """
    if u0_hat.ndim != 4 or u0_hat.shape[0] != 3:
        raise ValueError("u0_hat must have shape (3, N, N, N)")
    N = int(u0_hat.shape[1])
    if u0_hat.shape[2] != N or u0_hat.shape[3] != N:
        raise ValueError("u0_hat must be cubic in space")

    kx, ky, kz, k2 = kgrid_3d(N, L=2 * np.pi)
    k_mag = np.sqrt(k2)
    mask = dyadic_shell_mask(k_mag, int(j), k0=1.0).astype(float)

    out = np.array(u0_hat, dtype=complex, copy=True)
    out = out * mask[None, ...]
    out = project_div_free_3d(out, kx, ky, kz)
    if dealias:
        out = out * dealias_mask_3d(N)[None, ...]
        out = project_div_free_3d(out, kx, ky, kz)
    out[:, 0, 0, 0] = 0.0
    return out


def rescale_energy(u0_hat: np.ndarray, target_energy: float) -> np.ndarray:
    """Rescale a Fourier-space field to the requested physical energy."""
    if target_energy < 0:
        raise ValueError("target_energy must be nonnegative")
    if u0_hat.ndim != 4 or u0_hat.shape[0] != 3:
        raise ValueError("u0_hat must have shape (3, N, N, N)")
    current = _energy_from_u_hat(u0_hat)
    if current <= 1e-30 or target_energy == 0.0:
        return np.zeros_like(u0_hat)
    scale = np.sqrt(target_energy / current)
    return np.asarray(u0_hat, dtype=complex) * scale
