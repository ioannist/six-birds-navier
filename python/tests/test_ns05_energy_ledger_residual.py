import numpy as np

from nswave.ns2d_spectral import simulate_ns2d


def _cumulative_trapz(t: np.ndarray, y: np.ndarray) -> np.ndarray:
    out = np.zeros_like(y)
    if len(t) < 2:
        return out
    dt = t[1:] - t[:-1]
    out[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * dt)
    return out


def test_energy_ledger_residual_small():
    res = simulate_ns2d(
        N=32,
        dt=5e-4,
        t_max=0.05,
        nu=1e-3,
        mu=1e-4,
        alpha=2.0,
        forcing_hat_fn=None,
        dealias=True,
        record_every=1,
    )
    assert res.success is True

    t = res.t
    E = res.energy
    D = res.diss_nu + res.diss_mu

    assert np.all(np.isfinite(t))
    assert np.all(np.isfinite(E))
    assert np.all(np.isfinite(D))

    I = _cumulative_trapz(t, D)
    E0 = E[0]
    denom = max(E0, 1e-12)
    R = (E + I - E0) / denom

    resid_final = float(R[-1])
    max_abs_resid = float(np.max(np.abs(R)))

    assert abs(resid_final) <= 1e-2
    assert max_abs_resid <= 2e-2
