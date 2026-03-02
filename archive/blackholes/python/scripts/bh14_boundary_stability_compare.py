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

from bhwave.barrier import poeschl_teller_V  # noqa: E402
from bhwave.time_domain import ImpedanceBoundary, discrete_energy, simulate_wave_fdtd  # noqa: E402


def run_variant(boundary_update: str, boundary_bc: str) -> tuple[np.ndarray, np.ndarray, bool, int]:
    V0 = 1.0
    a = 1.0
    x_peak = 0.0

    x0 = -8.0
    x_max = 40.0
    x_obs = 20.0

    dx = 0.02
    dt = 0.018
    t_max = 120.0
    sponge_start = 30.0
    sponge_gamma_max = 1.0

    boundary = ImpedanceBoundary.from_debye(a0=0.2, a=[0.2], b=[1.0])

    x_src = 10.0
    sigma = 0.8

    def psi_init(x: np.ndarray) -> np.ndarray:
        return np.exp(-((x - x_src) ** 2) / (2.0 * sigma**2))

    def psi_t_init(x: np.ndarray) -> np.ndarray:
        base = psi_init(x)
        return -((x - x_src) / (sigma**2)) * base

    def V_func(x: np.ndarray) -> np.ndarray:
        return V0 / np.cosh(a * (x - x_peak)) ** 2

    record_full_every = 10
    result = simulate_wave_fdtd(
        V_func=V_func,
        x0=x0,
        x_max=x_max,
        dx=dx,
        dt=dt,
        t_max=t_max,
        boundary=boundary,
        x_obs=x_obs,
        sponge_start=sponge_start,
        sponge_gamma_max=sponge_gamma_max,
        psi_init=psi_init,
        psi_t_init=psi_t_init,
        record_full_every=record_full_every,
        boundary_update=boundary_update,
        boundary_bc=boundary_bc,
    )

    if result.psi_full is None or result.t_full is None or result.q_full is None:
        raise RuntimeError("full-field recording missing")

    V_grid = V_func(result.x)
    energies = []
    for i in range(1, len(result.t_full)):
        psi_curr = result.psi_full[i]
        psi_prev = result.psi_full[i - 1]
        q_snap = result.q_full[i]
        boundary_snap = ImpedanceBoundary(a0=boundary.a0, a=boundary.a, b=boundary.b, q=q_snap)
        dt_eff = dt * record_full_every
        energies.append(discrete_energy(psi_curr, psi_prev, V_grid, dx, dt_eff, boundary=boundary_snap))

    energies = np.asarray(energies, dtype=float)
    times = result.t_full[1:]
    finite_mask = np.isfinite(energies)
    has_nan = not np.all(finite_mask)
    finite_count = int(np.sum(finite_mask))
    return times, energies, has_nan, finite_count


def main() -> None:
    artifacts_dir = ROOT / "artifacts"
    artifacts_dir.mkdir(parents=True, exist_ok=True)

    times_euler, energy_euler, nan_euler, count_euler = run_variant("euler", "simple")
    times_exact, energy_exact, nan_exact, count_exact = run_variant("exact", "predictor_corrector")

    def max_rel_increase(energy: np.ndarray, times: np.ndarray, t_cut: float | None = None) -> float:
        if t_cut is None:
            mask = np.isfinite(energy)
        else:
            mask = (times >= t_cut) & np.isfinite(energy)
        finite = energy[mask]
        if finite.size == 0:
            return float("nan")
        E0 = float(finite[0])
        return float(np.max((finite - E0) / max(E0, 1e-12)))

    # Focus on drift after the packet has interacted with the inner boundary.
    t_cut = 23.0
    inc_euler_full = max_rel_increase(energy_euler, times_euler)
    inc_exact_full = max_rel_increase(energy_exact, times_exact)
    inc_euler = max_rel_increase(energy_euler, times_euler, t_cut=t_cut)
    inc_exact = max_rel_increase(energy_exact, times_exact, t_cut=t_cut)

    print("BH-14 boundary update stability compare")
    print(f"  max_rel_increase_euler={inc_euler:.6g}")
    print(f"  max_rel_increase_exact={inc_exact:.6g}")
    print(f"  max_rel_increase_euler_full={inc_euler_full:.6g}")
    print(f"  max_rel_increase_exact_full={inc_exact_full:.6g}")
    print(f"  max_rel_increase_window_start={t_cut:.1f}")
    print(f"  nan_euler={nan_euler}, nan_exact={nan_exact}")
    print(f"  finite_points_euler={count_euler}, finite_points_exact={count_exact}")

    plt.figure(figsize=(6.0, 4.0))
    plt.plot(times_euler, energy_euler, label="euler")
    plt.plot(times_exact, energy_exact, label="exact")
    plt.xlabel("t")
    plt.ylabel("Energy")
    plt.title("BH-14 boundary energy compare")
    plt.legend()
    plt.tight_layout()
    plot_path = artifacts_dir / "bh14_energy_compare.png"
    plt.savefig(plot_path)
    print(f"Saved plot: {plot_path}")


if __name__ == "__main__":
    main()
