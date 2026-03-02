import numpy as np

from nswave.ns2d_spectral import simulate_ns2d


def test_unforced_enstrophy_decreases():
    res = simulate_ns2d(
        N=32,
        dt=5e-4,
        t_max=0.05,
        nu=1e-3,
        mu=1e-4,
        alpha=2.0,
        forcing_hat_fn=None,
        record_every=10,
    )
    assert res.success is True

    Z = res.enstrophy
    assert np.all(np.isfinite(Z))
    max_increase = np.max(Z[1:] - Z[:-1])
    assert max_increase <= 1e-4 * Z[0] + 1e-12


def test_conservative_enstrophy_approx_constant():
    res = simulate_ns2d(
        N=32,
        dt=2e-4,
        t_max=0.01,
        nu=0.0,
        mu=0.0,
        alpha=2.0,
        forcing_hat_fn=None,
        record_every=5,
    )
    Z = res.enstrophy
    assert np.all(np.isfinite(Z))
    rel_err = abs(Z[-1] - Z[0]) / max(Z[0], 1e-12)
    assert rel_err <= 1e-2
