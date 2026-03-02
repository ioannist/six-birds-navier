import numpy as np

from nswave.ns3d_spectral import make_deterministic_u0_hat_3d, simulate_ns3d


def test_ns3d_dissipative_smoke():
    N = 16
    dt = 5e-4
    t_max = 0.01
    nu = 1e-3
    mu = 1e-4
    alpha = 2.0

    u0_hat = make_deterministic_u0_hat_3d(N, seed=0, k_init=6)
    res = simulate_ns3d(
        N=N,
        dt=dt,
        t_max=t_max,
        nu=nu,
        mu=mu,
        alpha=alpha,
        u0_hat=u0_hat,
        record_every=1,
    )

    assert res.success is True
    max_div = float(np.max(res.div_rms))
    assert max_div <= 1e-6
    assert res.energy[-1] <= res.energy[0] + 1e-10


def test_ns3d_conservative_short():
    N = 16
    dt = 2e-4
    t_max = 0.002
    nu = 0.0
    mu = 0.0
    alpha = 2.0

    u0_hat = make_deterministic_u0_hat_3d(N, seed=0, k_init=6)
    res = simulate_ns3d(
        N=N,
        dt=dt,
        t_max=t_max,
        nu=nu,
        mu=mu,
        alpha=alpha,
        u0_hat=u0_hat,
        record_every=1,
    )

    assert res.success is True
    rel_err = abs(res.energy[-1] - res.energy[0]) / max(res.energy[0], 1e-12)
    assert rel_err <= 5e-2
