#!/usr/bin/env python3
from __future__ import annotations

import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

from bhwave.conventions import R_from_Z, is_passive_R, is_passive_Z  # noqa: E402


def flux(omega: float, R: complex) -> float:
    return omega * (abs(R) ** 2 - 1.0)


def report_family(name: str, omega_grid: np.ndarray, Z_vals: np.ndarray) -> None:
    R_vals = np.array([R_from_Z(z) for z in Z_vals])
    min_re = float(np.min(Z_vals.real))
    max_absR = float(np.max(np.abs(R_vals)))
    max_flux = float(np.max([flux(omega, r) for omega, r in zip(omega_grid, R_vals)]))

    passive_min = min_re >= -1e-12
    pass_ok = (
        passive_min
        and max_absR <= 1.0 + 1e-9
        and max_flux <= 1e-9
    )
    status = "PASS" if pass_ok else "FAIL"

    print(f"{name} | min Re(Z)={min_re:.6g} | max |R|={max_absR:.6g} | max flux={max_flux:.6g} | {status}")


def identity_check(Z: complex) -> None:
    lhs = abs(1 + Z) ** 2 - abs(1 - Z) ** 2
    rhs = 4.0 * Z.real
    print(f"Z={Z}: |1+Z|^2 - |1-Z|^2 = {lhs:.6g}, 4 Re(Z) = {rhs:.6g}")


def main() -> None:
    omega_grid = np.logspace(-2, 2, 801)

    const_cases = [
        ("const Z=1.0", 1.0 + 0j),
        ("const Z=0.2+0.3j", 0.2 + 0.3j),
        ("const Z=0+0.5j", 0.0 + 0.5j),
        ("const Z=-0.2", -0.2 + 0j),
    ]

    print("Passivity table:")
    for name, Z in const_cases:
        Z_vals = np.full_like(omega_grid, Z, dtype=np.complex128)
        report_family(name, omega_grid, Z_vals)
        R = R_from_Z(Z)
        if name == "const Z=-0.2":
            print(f"  active example |R|={abs(R):.6g} (should exceed 1)")

    families = [
        ("rational a0=0.2,a1=0.7,b=3.0", 0.2, 0.7, 3.0),
        ("rational a0=0.0,a1=1.0,b=1.0", 0.0, 1.0, 1.0),
    ]
    for name, a0, a1, b in families:
        Z_vals = a0 + a1 / (1.0 - 1j * omega_grid / b)
        report_family(name, omega_grid, Z_vals)

    print("\nIdentity checks:")
    identity_check(0.2 + 0.3j)
    identity_check(0.0 + 0.5j)
    identity_check(-0.2 + 0j)

    # Quick internal consistency checks for mapping helpers.
    Z = 0.2 + 0.3j
    R = R_from_Z(Z)
    print("\nHelper checks:")
    print(f"is_passive_Z({Z}) -> {is_passive_Z(Z)}")
    print(f"is_passive_R({R}) -> {is_passive_R(R)}")
    print(f"flux(omega=2.0) -> {flux(2.0, R):.6g}")


if __name__ == "__main__":
    main()
