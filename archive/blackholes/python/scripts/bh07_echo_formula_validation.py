#!/usr/bin/env python3
from __future__ import annotations

import cmath
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

from bhwave.barrier import barrier_reflection_transmission, poeschl_teller_V  # noqa: E402
from bhwave.echo import cavity_H, predicted_R_out_x0  # noqa: E402
from bhwave.scattering import scattering_reflection  # noqa: E402


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
    R0 = 0.2 * cmath.exp(0.3j)

    omega_grid = np.linspace(0.5, 5.0, 101)

    print("Config:")
    print(f"  x0={x0}, x_ref={x_ref}, Dx={Dx}")
    print(f"  barrier: V0={V0}, a={a}, x_peak={x_peak}, x_min={x_min}, x_max={x_max}")
    print(f"  R0={R0}")
    print()

    abs_err_R = []
    rel_err_R = []
    rel_err_H = []
    h_count = 0

    sample_omegas = [0.5, 1.0, 2.0, 3.5, 5.0]

    for omega in omega_grid:
        barrier = barrier_reflection_transmission(
            float(omega),
            V,
            x_min=x_min,
            x_max=x_max,
            x_ref=x_ref,
        )
        r = barrier.R_barrier
        t = barrier.T_barrier

        full = scattering_reflection(
            float(omega),
            V,
            x0,
            x_max,
            R_boundary=R0,
        )
        R_full = full.R_out

        R_pred = predicted_R_out_x0(float(omega), R0, r, t, Dx)
        abs_err = abs(R_full - R_pred)
        abs_err_R.append(abs_err)
        rel_err_R.append(abs_err / max(1e-12, abs(R_full)))

        H_pred = cavity_H(float(omega), R0, r, Dx)
        if abs(t) >= 0.2:
            H_num = (R_full - r * cmath.exp(-2j * omega * Dx)) / (t * t * R0)
            rel_err = abs(H_num - H_pred) / max(1e-12, abs(H_pred))
            rel_err_H.append(rel_err)
            h_count += 1
        else:
            H_pred = None
            rel_err = None

    print("Summary:")
    print(f"  max_abs_err_R={max(abs_err_R):.6g}")
    print(f"  max_rel_err_R={max(rel_err_R):.6g}")
    if rel_err_H:
        print(f"  max_rel_err_H={max(rel_err_H):.6g}")
    else:
        print("  max_rel_err_H=nan")
    print(f"  H_points_used={h_count}")
    print()

    print("Sample table:")
    print("omega | R_full | R_pred | abs_err | |t| | rel_err_H")
    for omega in sample_omegas:
        barrier = barrier_reflection_transmission(
            float(omega),
            V,
            x_min=x_min,
            x_max=x_max,
            x_ref=x_ref,
        )
        r = barrier.R_barrier
        t = barrier.T_barrier
        full = scattering_reflection(
            float(omega),
            V,
            x0,
            x_max,
            R_boundary=R0,
        )
        R_full = full.R_out
        R_pred = predicted_R_out_x0(float(omega), R0, r, t, Dx)
        abs_err = abs(R_full - R_pred)
        rel_err_H = None
        if abs(t) >= 0.2:
            H_pred = cavity_H(float(omega), R0, r, Dx)
            H_num = (R_full - r * cmath.exp(-2j * omega * Dx)) / (t * t * R0)
            rel_err_H = abs(H_num - H_pred) / max(1e-12, abs(H_pred))
        rel_err_str = "skip" if rel_err_H is None else f"{rel_err_H:.6g}"
        print(
            f"{omega} | {R_full} | {R_pred} | {abs_err:.6g} | "
            f"{abs(t):.6g} | {rel_err_str}"
        )


if __name__ == "__main__":
    main()
