import numpy as np

from bhwave.time_domain import ImpedanceBoundary, simulate_wave_fdtd


def test_time_domain_smoke():
    def V_func(x: np.ndarray) -> np.ndarray:
        return np.zeros_like(x)

    x0 = 0.0
    x_max = 2.0
    dx = 0.02
    dt = 0.018
    t_max = 1.0
    boundary = ImpedanceBoundary.from_constant_Z(1.0)

    def psi_init(x: np.ndarray) -> np.ndarray:
        return np.exp(-((x - 1.0) ** 2) / (2.0 * 0.1**2))

    result = simulate_wave_fdtd(
        V_func=V_func,
        x0=x0,
        x_max=x_max,
        dx=dx,
        dt=dt,
        t_max=t_max,
        boundary=boundary,
        x_obs=1.0,
        psi_init=psi_init,
    )

    assert np.all(np.isfinite(result.psi_obs))
    assert len(result.t) == len(result.psi_obs)
