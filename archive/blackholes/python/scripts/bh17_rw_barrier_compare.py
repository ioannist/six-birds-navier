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

from bhwave.barrier import (  # noqa: E402
    barrier_reflection_transmission,
    fit_poeschl_teller_from_samples,
    regge_wheeler_V_rstar,
)


def main() -> None:
    M = 1.0
    ell = 2

    x_un = np.linspace(-80.0, 140.0, 22001)
    V_rw_un = regge_wheeler_V_rstar(x_un, M=M, ell=ell)

    i_peak = int(np.argmax(V_rw_un))
    x_peak = float(x_un[i_peak])
    x = x_un - x_peak

    x_peak_fit, V0_fit, a_fit = fit_poeschl_teller_from_samples(x, V_rw_un)
    V_pt_grid = V0_fit / np.cosh(a_fit * (x - x_peak_fit)) ** 2

    def Vrw_scalar(xq: float) -> float:
        return float(np.interp(xq, x, V_rw_un))

    def Vpt_scalar(xq: float) -> float:
        return float(np.interp(xq, x, V_pt_grid))

    x_min = float(x[0])
    x_max = float(x[-1])
    x_ref = 0.0

    omega_grid = np.linspace(0.2, 5.0, 61)

    R_rw = np.empty_like(omega_grid, dtype=complex)
    T_rw = np.empty_like(omega_grid, dtype=complex)
    R_pt = np.empty_like(omega_grid, dtype=complex)
    T_pt = np.empty_like(omega_grid, dtype=complex)

    for i, omega in enumerate(omega_grid):
        res_rw = barrier_reflection_transmission(
            float(omega),
            Vrw_scalar,
            x_min=x_min,
            x_max=x_max,
            x_ref=x_ref,
        )
        res_pt = barrier_reflection_transmission(
            float(omega),
            Vpt_scalar,
            x_min=x_min,
            x_max=x_max,
            x_ref=x_ref,
        )
        R_rw[i] = res_rw.R_barrier
        T_rw[i] = res_rw.T_barrier
        R_pt[i] = res_pt.R_barrier
        T_pt[i] = res_pt.T_barrier

    u_rw = np.abs(R_rw) ** 2 + np.abs(T_rw) ** 2 - 1.0
    u_pt = np.abs(R_pt) ** 2 + np.abs(T_pt) ** 2 - 1.0

    rms_R = float(np.sqrt(np.mean((np.abs(R_rw) - np.abs(R_pt)) ** 2)))

    print("BH-17 RW vs PT barrier compare")
    print(f"  RW peak x* (unshifted) = {x_peak:.6f}")
    print(f"  PT fit: V0={V0_fit:.6g}, a={a_fit:.6g}, x_peak_fit={x_peak_fit:.6g}")
    print(f"  max |unitarity residual| RW: {float(np.max(np.abs(u_rw))):.3e}")
    print(f"  max |unitarity residual| PT: {float(np.max(np.abs(u_pt))):.3e}")
    print(f"  RMS diff |R| (RW vs PT): {rms_R:.6g}")

    artifacts_dir = ROOT / "artifacts"
    artifacts_dir.mkdir(parents=True, exist_ok=True)
    plot_path = artifacts_dir / "bh17_rw_vs_pt.png"

    fig, axes = plt.subplots(2, 1, figsize=(6.4, 6.4), sharex=True)
    axes[0].plot(omega_grid, np.abs(R_rw), label="|R_RW|")
    axes[0].plot(omega_grid, np.abs(R_pt), label="|R_PT|", linestyle="--")
    axes[0].set_ylabel("|R|")
    axes[0].legend()

    axes[1].plot(omega_grid, np.abs(T_rw), label="|T_RW|")
    axes[1].plot(omega_grid, np.abs(T_pt), label="|T_PT|", linestyle="--")
    axes[1].set_xlabel("omega")
    axes[1].set_ylabel("|T|")
    axes[1].legend()

    fig.suptitle("BH-17: Regge-Wheeler vs Poeschl-Teller")
    fig.tight_layout()
    fig.savefig(plot_path)
    print(f"Saved plot: {plot_path}")


if __name__ == "__main__":
    main()
