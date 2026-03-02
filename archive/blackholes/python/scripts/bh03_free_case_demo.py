#!/usr/bin/env python3
from __future__ import annotations

import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

from bhwave.conventions import R_from_Z  # noqa: E402
from bhwave.scattering import scattering_reflection  # noqa: E402


def V0(x: float) -> float:
    return 0.0


def main() -> None:
    x0 = 0.0
    x_max = 5.0
    cases = [
        (0.7, 1.0 + 0j),
        (2.0, 0.2 + 0.3j),
        (4.5, 0.0 + 0.5j),
    ]

    print("omega | Z | R_expected | R_out | abs_error")
    for omega, Z in cases:
        res = scattering_reflection(omega, V0, x0, x_max, Z=Z)
        R_expected = R_from_Z(Z)
        err = abs(res.R_out - R_expected)
        print(f"{omega} | {Z} | {R_expected} | {res.R_out} | {err}")


if __name__ == "__main__":
    main()
