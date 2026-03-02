import numpy as np

from bhwave.time_domain import ImpedanceBoundary, discrete_energy, simulate_wave_fdtd


def _run(boundary_update: str):
    x0 = 0.0
    x_max = 2.0
    dx = 0.02
    dt = 0.018
    t_max = 2.0
    boundary = ImpedanceBoundary.from_debye(a0=0.2, a=[0.2], b=[1.0])

    def V_func(x: np.ndarray) -> np.ndarray:
        return np.zeros_like(x)

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
        record_full_every=1,
        boundary_update=boundary_update,
    )

    assert result.psi_full is not None
    assert result.t_full is not None
    assert result.q_full is not None

    V_grid = V_func(result.x)
    energies = []
    for i in range(1, len(result.t_full)):
        psi_curr = result.psi_full[i]
        psi_prev = result.psi_full[i - 1]
        q_snap = result.q_full[i]
        boundary_snap = ImpedanceBoundary(a0=boundary.a0, a=boundary.a, b=boundary.b, q=q_snap)
        energies.append(discrete_energy(psi_curr, psi_prev, V_grid, dx, dt, boundary=boundary_snap))

    energies = np.asarray(energies, dtype=float)
    assert np.all(np.isfinite(energies))
    E0 = float(energies[0])
    max_rel_increase = float(np.max((energies - E0) / max(E0, 1e-12)))
    return max_rel_increase


def test_exact_update_not_worse():
    inc_euler = _run("euler")
    inc_exact = _run("exact")
    assert inc_exact <= inc_euler + 1e-6
