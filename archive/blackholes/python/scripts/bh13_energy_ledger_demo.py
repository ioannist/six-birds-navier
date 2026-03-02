#!/usr/bin/env python3
from __future__ import annotations

import pathlib
import sys

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

from bhwave.time_domain import ImpedanceBoundary, discrete_energy, simulate_wave_fdtd  # noqa: E402


def main() -> None:
    artifacts_dir = ROOT / "artifacts"
    artifacts_dir.mkdir(parents=True, exist_ok=True)

    x0 = 0.0
    x_max = 10.0
    dx = 0.02
    dt = 0.012
    t_max = 50.0
    x_obs = 5.0

    boundary = ImpedanceBoundary.from_debye(a0=0.2, a=[0.2], b=[1.0])

    x_src = 6.0
    sigma = 0.5

    def V_func(x: np.ndarray) -> np.ndarray:
        return np.zeros_like(x)

    def psi_init(x: np.ndarray) -> np.ndarray:
        return np.zeros_like(x)

    def psi_t_init(x: np.ndarray) -> np.ndarray:
        return np.exp(-((x - x_src) ** 2) / (2.0 * sigma**2))

    result = simulate_wave_fdtd(
        V_func=V_func,
        x0=x0,
        x_max=x_max,
        dx=dx,
        dt=dt,
        t_max=t_max,
        boundary=boundary,
        x_obs=x_obs,
        psi_init=psi_init,
        psi_t_init=psi_t_init,
        record_full_every=1,
    )

    if result.psi_full is None or result.t_full is None or result.q_full is None:
        raise RuntimeError("full-field recording missing")

    V_grid = V_func(result.x)
    energies = []
    times = []
    for i in range(1, len(result.t_full)):
        psi_curr = result.psi_full[i]
        psi_prev = result.psi_full[i - 1]
        q_snap = result.q_full[i]
        boundary_snap = ImpedanceBoundary(a0=boundary.a0, a=boundary.a, b=boundary.b, q=q_snap)
        E = discrete_energy(
            psi_curr,
            psi_prev,
            V_grid,
            dx,
            dt,
            boundary=boundary_snap,
        )
        energies.append(E)
        times.append(result.t_full[i])

    energies = np.asarray(energies, dtype=float)
    times = np.asarray(times, dtype=float)
    E0 = float(energies[0])
    Emin = float(np.min(energies))
    Efinal = float(energies[-1])
    max_rel_increase = float(np.max((energies - E0) / max(E0, 1e-12)))

    print("Config:")
    print(f"  dx={dx}, dt={dt}, CFL={dt/dx:.3f}")
    print(f"  boundary: a0={boundary.a0}, a={boundary.a.tolist()}, b={boundary.b.tolist()}")
    print()
    print(f"E0={E0:.6g}")
    print(f"Emin={Emin:.6g}")
    print(f"Efinal={Efinal:.6g}")
    print(f"max_rel_increase={max_rel_increase:.6g}")

    plt.figure(figsize=(6.0, 4.0))
    plt.plot(times, energies, label="E(t)")
    plt.xlabel("t")
    plt.ylabel("Energy")
    plt.title("BH-13 discrete energy ledger")
    plt.tight_layout()
    plot_path = artifacts_dir / "bh13_energy.png"
    plt.savefig(plot_path)
    print(f"Saved plot: {plot_path}")


if __name__ == "__main__":
    main()
