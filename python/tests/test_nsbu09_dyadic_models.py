import numpy as np

from nswave.dyadic_models import simulate_kp


def test_inviscid_energy_conservation() -> None:
    J = 20
    dt = 1e-4
    t_max = 0.01
    u0 = np.zeros(J)
    u0[0] = 1.0
    u0[1] = 0.5
    res = simulate_kp(J=J, dt=dt, t_max=t_max, nu=0.0, mu=0.0, alpha=1.25, u0=u0)
    E = res["E_hist"]
    rel_drift = abs(E[-1] - E[0]) / max(E[0], 1e-12)
    assert rel_drift <= 1e-3


def test_passive_damping_decreases_energy() -> None:
    J = 20
    dt = 1e-4
    t_max = 0.01
    u0 = np.zeros(J)
    u0[0] = 1.0
    u0[1] = 0.5
    res = simulate_kp(J=J, dt=dt, t_max=t_max, nu=0.0, mu=1e-4, alpha=1.25, u0=u0)
    E = res["E_hist"]
    assert E[-1] <= E[0] + 1e-6


def test_active_sgs_not_strictly_decreasing() -> None:
    J = 20
    dt = 1e-4
    t_max = 0.01
    u0 = np.zeros(J)
    u0[0] = 1.0
    u0[1] = 0.5
    res = simulate_kp(J=J, dt=dt, t_max=t_max, nu=0.0, mu=-1e-4, alpha=1.25, u0=u0)
    E = res["E_hist"]
    assert E[-1] >= E[0] - 1e-6
