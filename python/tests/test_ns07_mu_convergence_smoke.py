import numpy as np

from nswave.ns2d_spectral import make_deterministic_omega0_hat, simulate_ns2d
from nswave.operators import kgrid_2d


def _low_k_error(omega_hat_base, omega_hat_mu, k2, L, k_low):
    k_mag = np.sqrt(k2)
    mask = k_mag <= k_low
    omega_low_base = np.fft.ifft2(omega_hat_base * mask).real
    omega_low_mu = np.fft.ifft2(omega_hat_mu * mask).real
    norm_base = np.sqrt(np.mean(omega_low_base**2) * L * L)
    norm_err = np.sqrt(np.mean((omega_low_mu - omega_low_base) ** 2) * L * L)
    return float(norm_err / max(norm_base, 1e-12))


def test_mu_convergence_ordering():
    N = 32
    L = 2 * np.pi
    dt = 5e-4
    t_max = 0.1
    nu = 1e-3
    alpha = 2.0
    dealias = True
    k_low = 4

    omega0_hat = make_deterministic_omega0_hat(N, dealias=dealias)
    _, _, k2 = kgrid_2d(N, L=L)

    base = simulate_ns2d(
        N=N,
        dt=dt,
        t_max=t_max,
        nu=nu,
        mu=0.0,
        alpha=alpha,
        omega0_hat=omega0_hat,
        forcing_hat_fn=None,
        L=L,
        dealias=dealias,
        record_every=10,
    )

    hi = simulate_ns2d(
        N=N,
        dt=dt,
        t_max=t_max,
        nu=nu,
        mu=1e-2,
        alpha=alpha,
        omega0_hat=omega0_hat,
        forcing_hat_fn=None,
        L=L,
        dealias=dealias,
        record_every=10,
    )

    lo = simulate_ns2d(
        N=N,
        dt=dt,
        t_max=t_max,
        nu=nu,
        mu=1e-3,
        alpha=alpha,
        omega0_hat=omega0_hat,
        forcing_hat_fn=None,
        L=L,
        dealias=dealias,
        record_every=10,
    )

    err_hi = _low_k_error(base.omega_hat_final, hi.omega_hat_final, k2, L, k_low)
    err_lo = _low_k_error(base.omega_hat_final, lo.omega_hat_final, k2, L, k_low)

    assert np.isfinite(err_hi)
    assert np.isfinite(err_lo)
    assert err_lo <= err_hi + 1e-6
