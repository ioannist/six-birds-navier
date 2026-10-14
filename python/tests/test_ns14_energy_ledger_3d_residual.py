import numpy as np

from nswave.ns3d_spectral import make_deterministic_u0_hat_3d, simulate_ns3d


def _cumulative_trapz(t: np.ndarray, y: np.ndarray) -> np.ndarray:
    out = np.zeros_like(y)
    if len(t) < 2:
        return out
    dt = t[1:] - t[:-1]
    out[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * dt)
    return out


def test_ns14_energy_ledger_residual():
    N = 16
    L = 2 * np.pi
    dt = 5e-4
    t_max = 0.01
    nu = 1e-3
    mu = 1e-4
    alpha = 2.0
    seed = 0
    k_init = 6

    u0_hat = make_deterministic_u0_hat_3d(
        N,
        L=L,
        seed=seed,
        k_init=k_init,
        dealias=True,
    )

    res = simulate_ns3d(
        N=N,
        dt=dt,
        t_max=t_max,
        nu=nu,
        mu=mu,
        alpha=alpha,
        u0_hat=u0_hat,
        forcing_hat_fn=None,
        L=L,
        dealias=True,
        record_every=1,
    )

    assert res.success is True
    assert np.all(np.isfinite(res.energy))
    assert np.all(np.isfinite(res.diss_nu))
    assert np.all(np.isfinite(res.diss_mu))

    t = res.t
    E = res.energy
    D = res.diss_nu + res.diss_mu
    I = _cumulative_trapz(t, D)

    E0 = E[0]
    denom = max(E0, 1e-12)
    R = (E + I - E0) / denom

    resid_final = float(R[-1])
    max_abs_resid = float(np.max(np.abs(R)))

    assert abs(resid_final) <= 1e-2
    assert max_abs_resid <= 2e-2
    assert E[-1] <= E[0] + 1e-10
