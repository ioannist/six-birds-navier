import numpy as np

from nswave.operators import kgrid_2d


def _vorticity_to_velocity_hat(omega_hat, kx, ky, k2):
    psi_hat = np.zeros_like(omega_hat)
    mask = k2 != 0
    psi_hat[mask] = omega_hat[mask] / k2[mask]
    u_hat = np.empty((2,) + omega_hat.shape, dtype=complex)
    u_hat[0] = 1j * ky * psi_hat
    u_hat[1] = -1j * kx * psi_hat
    return u_hat


def _rate(u_hat, c_hat, L):
    u0 = np.fft.ifft2(u_hat[0]).real
    u1 = np.fft.ifft2(u_hat[1]).real
    c0 = np.fft.ifft2(c_hat[0]).real
    c1 = np.fft.ifft2(c_hat[1]).real
    return float(np.mean(u0 * c0 + u1 * c1) * L * L)


def test_good_vs_anti_diffusion_rates():
    L = 2 * np.pi
    N = 16
    M = 5
    nu0 = 1e-2

    kx, ky, k2 = kgrid_2d(N, L=L)
    rng = np.random.default_rng(0)

    for _ in range(M):
        omega = rng.standard_normal((N, N))
        omega_hat = np.fft.fft2(omega)
        u_hat = _vorticity_to_velocity_hat(omega_hat, kx, ky, k2)

        c_good = -nu0 * k2 * u_hat
        c_good[:, 0, 0] = 0.0
        c_anti = nu0 * k2 * u_hat
        c_anti[:, 0, 0] = 0.0

        rate_good = _rate(u_hat, c_good, L)
        rate_anti = _rate(u_hat, c_anti, L)

        assert rate_good <= 1e-10
        assert rate_anti >= -1e-10
