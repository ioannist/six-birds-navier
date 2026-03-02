#!/usr/bin/env python3
from __future__ import annotations

import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

from bhwave.conventions import (  # noqa: E402
    R_from_Z,
    Z_from_R,
    Y_from_Z,
    R_from_Y,
    Z_from_Y,
    Y_from_R,
)


def main() -> None:
    omega = 2.0

    print("Formulas:")
    print("  R(Z) = (1 - Z) / (1 + Z)")
    print("  Z(R) = (1 - R) / (1 + R)")
    print("  Y(Z) = -i * omega * Z")
    print("  R(Y) = (i * omega + Y) / (i * omega - Y)")
    print("  Z(Y) = i * Y / omega")
    print("  Y(R) = -i * omega * (1 - R) / (1 + R)")
    print(f"  omega = {omega}")
    print()

    checks = []

    Z = 1 + 0j
    R = R_from_Z(Z)
    checks.append(("Z=1 -> R=0", R))

    Z = 0 + 0j
    R = R_from_Z(Z)
    checks.append(("Z=0 -> R=1", R))

    Z = 1j * 0.5
    R = R_from_Z(Z)
    checks.append(("Z=0+0.5j -> |R|=1", abs(R)))

    Z = 0.2 + 0.3j
    R = R_from_Z(Z)
    checks.append(("Z=0.2+0.3j -> |R|<1, Re(Z)>0", (abs(R), Z.real)))

    Z = -0.2 + 0j
    R = R_from_Z(Z)
    checks.append(("Z=-0.2 -> |R|>1", abs(R)))

    R = 0 + 0j
    Y = Y_from_R(R, omega)
    checks.append(("R=0 -> Y=-i*omega", Y))

    R = 1 + 0j
    Y = Y_from_R(R, omega)
    checks.append(("R=1 -> Y=0", Y))

    Z = 1 + 0j
    R = R_from_Z(Z)
    Y = Y_from_Z(Z, omega)
    checks.append(("Z=1 -> R=0, Y=-i*omega", (R, Y)))

    print("Sanity checks:")
    for label, value in checks:
        print(f"  {label}: {value}")

    print()
    samples = [0.2 + 0.1j, 1.0 + 0.0j, 2.0 + 1.0j]
    max_err = 0.0
    for Z in samples:
        Z2 = Z_from_R(R_from_Z(Z))
        err = abs(Z2 - Z)
        if err > max_err:
            max_err = err
    print(f"Round-trip Z -> R -> Z max error: {max_err}")

    print()
    Z = 1 + 0j
    Y = Y_from_Z(Z, omega)
    Z_back = Z_from_Y(Y, omega)
    print(f"Y(Z=1) = {Y}, Z(Y) = {Z_back}")

    print()
    R_samples = [0 + 0j, 0.3 + 0.2j, 0.6 - 0.8j, 0.9 + 0j, -0.5 + 0j]
    max_err_R = 0.0
    for R in R_samples:
        R2 = R_from_Y(Y_from_R(R, omega), omega)
        err = abs(R2 - R)
        if err > max_err_R:
            max_err_R = err
    print(f"Round-trip R -> Y -> R max error: {max_err_R}")


if __name__ == "__main__":
    main()
