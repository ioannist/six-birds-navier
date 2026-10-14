#!/usr/bin/env python3
from __future__ import annotations

import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

from bhwave.barrier import barrier_reflection_transmission, poeschl_teller_V  # noqa: E402


def main() -> None:
    V0 = 1.0
    a = 1.0
    x_peak = 0.0

    def V(x: float) -> float:
        return poeschl_teller_V(x, V0, a, x_peak)

    x_min = -12.0
    x_max = 12.0
    x_ref = 0.0
    omegas = [0.3, 0.5, 1.0, 2.0, 5.0]

    print("omega | |R_barrier| | |T_barrier| | |R|^2+|T|^2-1")
    for omega in omegas:
        res = barrier_reflection_transmission(
            omega,
            V,
            x_min=x_min,
            x_max=x_max,
            x_ref=x_ref,
        )
        mag_R = abs(res.R_barrier)
        mag_T = abs(res.T_barrier)
        residual = mag_R**2 + mag_T**2 - 1.0
        print(f"{omega} | {mag_R:.6g} | {mag_T:.6g} | {residual:.6g}")


if __name__ == "__main__":
    main()
