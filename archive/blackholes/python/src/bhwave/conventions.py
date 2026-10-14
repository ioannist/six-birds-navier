"""BH-02 conventions used by the boundary mapping helpers.

Time dependence is exp(-i omega t) with y = x - x0, boundary at y = 0.
Near the boundary: psi = A_in e^{-i omega y} + A_out e^{+i omega y}.
Define R = A_out / A_in and psi'(0) = -i omega Z psi(0).
Then Y = psi'/psi = -i omega Z, with R = (1 - Z)/(1 + Z).
"""

from __future__ import annotations


def R_from_Z(Z: complex) -> complex:
    if 1 + Z == 0:
        raise ValueError("R_from_Z undefined for Z = -1")
    return (1 - Z) / (1 + Z)


def Z_from_R(R: complex) -> complex:
    if 1 + R == 0:
        raise ValueError("Z_from_R undefined for R = -1")
    return (1 - R) / (1 + R)


def Y_from_Z(Z: complex, omega: float) -> complex:
    return -1j * omega * Z


def R_from_Y(Y: complex, omega: float) -> complex:
    if omega == 0:
        raise ValueError("omega must be nonzero")
    if 1j * omega - Y == 0:
        raise ValueError("R_from_Y undefined for Y = i*omega")
    return (1j * omega + Y) / (1j * omega - Y)


def Z_from_Y(Y: complex, omega: float) -> complex:
    if omega == 0:
        raise ValueError("omega must be nonzero")
    return 1j * Y / omega


def Y_from_R(R: complex, omega: float) -> complex:
    if omega == 0:
        raise ValueError("omega must be nonzero")
    if 1 + R == 0:
        raise ValueError("Y_from_R undefined for R = -1")
    return -1j * omega * (1 - R) / (1 + R)


def is_passive_Z(Z: complex, tol: float = 1e-12) -> bool:
    return Z.real >= -tol


def is_passive_R(R: complex, tol: float = 1e-12) -> bool:
    return abs(R) <= 1 + tol
