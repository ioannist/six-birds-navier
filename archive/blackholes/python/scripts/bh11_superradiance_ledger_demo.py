#!/usr/bin/env python3
from __future__ import annotations

import cmath
import pathlib
import sys

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

from bhwave.barrier import barrier_reflection_transmission, poeschl_teller_V  # noqa: E402
from bhwave.superradiance import growth_rate_from_gain, pole_estimate_from_G, round_trip_gain  # noqa: E402


def main() -> None:
    artifacts_dir = ROOT / "artifacts"
    artifacts_dir.mkdir(parents=True, exist_ok=True)

    V0 = 1.0
    a = 1.0
    x_peak = 0.0
    x_ref = 0.0
    x_min = -12.0
    x_max = 12.0

    x0 = -8.0
    Dx = x_ref - x0

    alpha = 0.25
    omega_s = 0.5
    sigma = 0.2
    theta0 = 0.0

    def V(x: float) -> float:
        return poeschl_teller_V(x, V0, a, x_peak)

    omega_grid = np.linspace(0.2, 5.0, 241)
    r_vals = np.zeros_like(omega_grid, dtype=complex)
    gain_vals = np.zeros_like(omega_grid, dtype=float)

    for i, omega in enumerate(omega_grid):
        res = barrier_reflection_transmission(
            float(omega),
            V,
            x_min=x_min,
            x_max=x_max,
            x_ref=x_ref,
        )
        r = res.R_barrier
        rho = 1.0 + alpha * np.exp(-0.5 * ((omega - omega_s) / sigma) ** 2)
        R0 = rho * np.exp(1j * theta0)
        r_vals[i] = r
        gain_vals[i] = round_trip_gain(R0, r)

    idx = int(np.argmax(gain_vals))
    omega_star = float(omega_grid[idx])
    rho_star = 1.0 + alpha * np.exp(-0.5 * ((omega_star - omega_s) / sigma) ** 2)
    R0_star = rho_star * np.exp(1j * theta0)
    G_star = R0_star * r_vals[idx]
    g_star = abs(G_star)

    growth_rate = growth_rate_from_gain(g_star, Dx)
    n = int(round((2.0 * Dx * omega_star + cmath.phase(G_star)) / (2.0 * np.pi)))
    omega_pole = pole_estimate_from_G(G_star, Dx, n)

    print("Config:")
    print(f"  Dx={Dx}")
    print(f"  barrier: V0={V0}, a={a}, x_peak={x_peak}, x_min={x_min}, x_max={x_max}")
    print(f"  active band: alpha={alpha}, omega_s={omega_s}, sigma={sigma}, theta0={theta0}")
    print()

    print(f"max gain g*={g_star:.6g} at omega*={omega_star:.6g}")
    print(f"growth_rate={growth_rate:.6g}")
    print(f"omega_pole_est={omega_pole}")

    amps = [1.0]
    for _ in range(9):
        amps.append(amps[-1] * g_star)
    print("Echo amplitude sequence A0..A9:")
    print("  ", [float(a) for a in amps])

    plt.figure(figsize=(6.0, 4.0))
    plt.plot(omega_grid, gain_vals, label="gain |R0*r|")
    plt.axhline(1.0, color="red", linestyle="--", linewidth=1)
    plt.xlabel("omega")
    plt.ylabel("gain")
    plt.title("BH-11 toy superradiant gain")
    plt.tight_layout()
    plot_path = artifacts_dir / "bh11_superradiance_gain.png"
    plt.savefig(plot_path)
    print(f"Saved plot: {plot_path}")


if __name__ == "__main__":
    main()
