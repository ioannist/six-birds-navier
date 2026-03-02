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
from bhwave.echo import cavity_H, predicted_R_out_x0  # noqa: E402
from bhwave.scattering import scattering_reflection  # noqa: E402


def Z_debye(omega: float, a0: float, a1: float, b1: float) -> complex:
    return a0 + a1 / (1.0 - 1j * omega / b1)


def compute_point(
    omega: float,
    *,
    V,
    x0: float,
    x_min: float,
    x_max: float,
    x_ref: float,
    Dx: float,
    a0: float,
    a1: float,
    b1: float,
) -> dict[str, complex | float | None]:
    Z = Z_debye(omega, a0, a1, b1)
    R0 = (1.0 - Z) / (1.0 + Z)

    barrier = barrier_reflection_transmission(
        omega,
        V,
        x_min=x_min,
        x_max=x_max,
        x_ref=x_ref,
    )
    r = barrier.R_barrier
    t = barrier.T_barrier

    full = scattering_reflection(
        omega,
        V,
        x0,
        x_max,
        Z=Z,
    )
    R_full = full.R_out

    R_pred = predicted_R_out_x0(omega, R0, r, t, Dx)
    abs_err = abs(R_full - R_pred)
    rel_err = abs_err / max(1e-12, abs(R_full))

    rel_err_H = None
    if abs(t) >= 0.2 and abs(R0) >= 0.05:
        H_pred = cavity_H(omega, R0, r, Dx)
        H_num = (R_full - r * cmath.exp(-2j * omega * Dx)) / (t * t * R0)
        rel_err_H = abs(H_num - H_pred) / max(1e-12, abs(H_pred))

    return {
        "omega": omega,
        "Z": Z,
        "R0": R0,
        "r": r,
        "t": t,
        "R_full": R_full,
        "R_pred": R_pred,
        "abs_err": abs_err,
        "rel_err": rel_err,
        "rel_err_H": rel_err_H,
    }


def run_case(
    name: str,
    *,
    a0: float,
    a1: float,
    b1: float,
    omega_grid: np.ndarray,
    V,
    x0: float,
    x_min: float,
    x_max: float,
    x_ref: float,
    Dx: float,
) -> tuple[dict[str, float], list[float]]:
    max_abs_err_R = 0.0
    max_rel_err_R = 0.0
    max_rel_err_H = 0.0
    h_points = 0
    abs_err_series: list[float] = []

    for omega in omega_grid:
        point = compute_point(
            float(omega),
            V=V,
            x0=x0,
            x_min=x_min,
            x_max=x_max,
            x_ref=x_ref,
            Dx=Dx,
            a0=a0,
            a1=a1,
            b1=b1,
        )
        abs_err = float(point["abs_err"])
        rel_err = float(point["rel_err"])
        abs_err_series.append(abs_err)
        max_abs_err_R = max(max_abs_err_R, abs_err)
        max_rel_err_R = max(max_rel_err_R, rel_err)
        if point["rel_err_H"] is not None:
            rel_err_H = float(point["rel_err_H"])
            max_rel_err_H = max(max_rel_err_H, rel_err_H)
            h_points += 1

    metrics = {
        "max_abs_err_R": max_abs_err_R,
        "max_rel_err_R": max_rel_err_R,
        "max_rel_err_H": max_rel_err_H,
        "H_points_used": float(h_points),
    }
    print(f"Case {name}:")
    print(f"  max_abs_err_R={max_abs_err_R:.6g}")
    print(f"  max_rel_err_R={max_rel_err_R:.6g}")
    print(f"  max_rel_err_H={max_rel_err_H:.6g}")
    print(f"  H_points_used={h_points}")
    print()

    return metrics, abs_err_series


def main() -> None:
    V0 = 1.0
    a = 1.0
    x_ref = 0.0
    x_peak = 0.0

    def V(x: float) -> float:
        return poeschl_teller_V(x, V0, a, x_peak)

    x0 = -8.0
    x_min = -12.0
    x_max = 12.0
    Dx = x_ref - x0

    omega_grid = np.linspace(0.5, 5.0, 101)

    caseA = {"name": "A", "a0": 0.2, "a1": 0.5, "b1": 2.0}
    caseB = {"name": "B", "a0": 0.1, "a1": 0.8, "b1": 0.7}

    print("Config:")
    print(f"  x0={x0}, x_ref={x_ref}, Dx={Dx}")
    print(f"  barrier: V0={V0}, a={a}, x_peak={x_peak}, x_min={x_min}, x_max={x_max}")
    print()

    metrics_A, abs_err_A = run_case(
        caseA["name"],
        a0=caseA["a0"],
        a1=caseA["a1"],
        b1=caseA["b1"],
        omega_grid=omega_grid,
        V=V,
        x0=x0,
        x_min=x_min,
        x_max=x_max,
        x_ref=x_ref,
        Dx=Dx,
    )
    run_case(
        caseB["name"],
        a0=caseB["a0"],
        a1=caseB["a1"],
        b1=caseB["b1"],
        omega_grid=omega_grid,
        V=V,
        x0=x0,
        x_min=x_min,
        x_max=x_max,
        x_ref=x_ref,
        Dx=Dx,
    )

    print("Sample table (case A):")
    print("omega | Z | R0 | R_full | R_pred | abs_err | |t| | rel_err_H")
    for omega in [0.5, 1.0, 2.0, 3.5, 5.0]:
        point = compute_point(
            float(omega),
            V=V,
            x0=x0,
            x_min=x_min,
            x_max=x_max,
            x_ref=x_ref,
            Dx=Dx,
            a0=caseA["a0"],
            a1=caseA["a1"],
            b1=caseA["b1"],
        )
        rel_err_H = point["rel_err_H"]
        rel_err_str = "skip" if rel_err_H is None else f"{rel_err_H:.6g}"
        print(
            f"{omega} | {point['Z']} | {point['R0']} | {point['R_full']} | {point['R_pred']} | "
            f"{point['abs_err']:.6g} | {abs(point['t']):.6g} | {rel_err_str}"
        )

    artifacts_dir = ROOT / "artifacts"
    artifacts_dir.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(6.0, 4.0))
    plt.plot(omega_grid, abs_err_A, label="|R_full - R_pred|")
    plt.xlabel("omega")
    plt.ylabel("abs error")
    plt.title("BH-15 frequency-dependent echo error (case A)")
    plt.tight_layout()
    plot_path = artifacts_dir / "bh15_freqdep_echo_error.png"
    plt.savefig(plot_path)
    print(f"Saved plot: {plot_path}")


if __name__ == "__main__":
    main()
