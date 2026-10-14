import numpy as np

from nswave.dyadic import frontier_scale, shell_energies
from nswave.operators import kgrid_3d


def test_parseval_shell_energy_sum() -> None:
    N = 16
    L = 2 * np.pi
    rng = np.random.default_rng(0)
    u = rng.standard_normal((3, N, N, N))
    u_hat = np.zeros_like(u, dtype=complex)
    for i in range(3):
        u_hat[i] = np.fft.fftn(u[i])
    u_hat[:, 0, 0, 0] = 0.0
    u_phys = [np.fft.ifftn(u_hat[i]).real for i in range(3)]

    energy_phys = 0.5 * np.mean(u_phys[0] ** 2 + u_phys[1] ** 2 + u_phys[2] ** 2) * L**3

    kx, ky, kz, k2 = kgrid_3d(N, L=L)
    k_mag = np.sqrt(k2)
    E_shell = shell_energies(u_hat, k_mag, L=L)
    energy_shell = float(np.sum(E_shell))

    assert np.isclose(energy_shell, energy_phys, rtol=1e-10, atol=1e-12)


def test_frontier_scale_sanity() -> None:
    E = np.array([0.1, 0.05, 0.01])
    assert frontier_scale(E, 0.02) == 1
    assert frontier_scale(E, 0.2) == -1
