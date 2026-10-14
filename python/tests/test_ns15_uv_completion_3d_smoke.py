import numpy as np
import pytest

from nswave.ns3d_spectral import make_deterministic_u0_hat_3d, simulate_ns3d
from nswave.operators import kgrid_3d


@pytest.mark.slow
def test_ns15_uv_completion_tail_ratio():
    N = 16
    L = 2 * np.pi
    dt = 5e-4
    t_max = 0.02
    nu = 1e-3
    alpha = 2.0
    mu_base = 0.0
    mu_comp = 1e-2
    k_low = 3
    k0 = 5

    u0_hat = make_deterministic_u0_hat_3d(N, L=L, seed=0, k_init=6, dealias=True)

    res0 = simulate_ns3d(
        N=N,
        dt=dt,
        t_max=t_max,
        nu=nu,
        mu=mu_base,
        alpha=alpha,
        u0_hat=u0_hat,
        dealias=True,
        record_every=1,
    )
    res1 = simulate_ns3d(
        N=N,
        dt=dt,
        t_max=t_max,
        nu=nu,
        mu=mu_comp,
        alpha=alpha,
        u0_hat=u0_hat,
        dealias=True,
        record_every=1,
    )

    kx, ky, kz, k2 = kgrid_3d(N, L=L)
    k_mag = np.sqrt(k2)
    shell = np.floor(k_mag + 0.5).astype(int)
    e0 = 0.5 * (
        np.abs(res0.u_hat_final[0]) ** 2
        + np.abs(res0.u_hat_final[1]) ** 2
        + np.abs(res0.u_hat_final[2]) ** 2
    )
    e1 = 0.5 * (
        np.abs(res1.u_hat_final[0]) ** 2
        + np.abs(res1.u_hat_final[1]) ** 2
        + np.abs(res1.u_hat_final[2]) ** 2
    )
    kmax = int(shell.max())
    e0_shell = np.bincount(shell.ravel(), weights=e0.ravel(), minlength=kmax + 1)
    e1_shell = np.bincount(shell.ravel(), weights=e1.ravel(), minlength=kmax + 1)

    low0 = float(np.sum(e0_shell[: k_low + 1]))
    low1 = float(np.sum(e1_shell[: k_low + 1]))
    high0 = float(np.sum(e0_shell[k0:]))
    high1 = float(np.sum(e1_shell[k0:]))

    ratio_low = low1 / max(low0, 1e-30)
    ratio_high = high1 / max(high0, 1e-30)

    assert ratio_low > 0.5
    assert ratio_high < 0.9
