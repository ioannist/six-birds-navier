import numpy as np
import pytest

from nswave.ns3d_spectral import make_deterministic_u0_hat_3d, simulate_ns3d
from nswave.operators import kgrid_3d


@pytest.mark.slow
def test_ns16_mu_convergence_3d_smoke():
    N = 16
    L = 2 * np.pi
    dt = 5e-4
    t_max = 0.02
    nu = 1e-3
    alpha = 2.0
    mu_hi = 1e-2
    mu_lo = 1e-3
    k_init = 6
    k_low = 3

    u0_hat = make_deterministic_u0_hat_3d(N, L=L, seed=0, k_init=k_init, dealias=True)
    _, _, _, k2 = kgrid_3d(N, L=L)

    base = simulate_ns3d(
        N=N,
        dt=dt,
        t_max=t_max,
        nu=nu,
        mu=0.0,
        alpha=alpha,
        u0_hat=u0_hat,
        dealias=True,
        record_every=1,
    )

    hi = simulate_ns3d(
        N=N,
        dt=dt,
        t_max=t_max,
        nu=nu,
        mu=mu_hi,
        alpha=alpha,
        u0_hat=u0_hat,
        dealias=True,
        record_every=1,
    )

    lo = simulate_ns3d(
        N=N,
        dt=dt,
        t_max=t_max,
        nu=nu,
        mu=mu_lo,
        alpha=alpha,
        u0_hat=u0_hat,
        dealias=True,
        record_every=1,
    )

    k_mag = np.sqrt(k2)
    mask = k_mag <= k_low

    def low_k_err(u_hat):
        u = [np.fft.ifftn(u_hat[i] * mask).real for i in range(3)]
        norm = np.sqrt(np.mean(u[0] ** 2 + u[1] ** 2 + u[2] ** 2) * L**3)
        return u, norm

    u_base, norm_base = low_k_err(base.u_hat_final)
    u_hi, _ = low_k_err(hi.u_hat_final)
    u_lo, _ = low_k_err(lo.u_hat_final)

    diff_hi = (u_hi[0] - u_base[0]) ** 2 + (u_hi[1] - u_base[1]) ** 2 + (u_hi[2] - u_base[2]) ** 2
    diff_lo = (u_lo[0] - u_base[0]) ** 2 + (u_lo[1] - u_base[1]) ** 2 + (u_lo[2] - u_base[2]) ** 2

    err_hi = np.sqrt(np.mean(diff_hi) * L**3) / max(norm_base, 1e-12)
    err_lo = np.sqrt(np.mean(diff_lo) * L**3) / max(norm_base, 1e-12)

    assert np.isfinite(err_hi)
    assert np.isfinite(err_lo)
    assert err_lo <= err_hi + 1e-10
