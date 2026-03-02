"""Frequency-domain 1D scattering solver with inner boundary response."""

from __future__ import annotations

from dataclasses import dataclass
import cmath
from typing import Callable

from scipy.integrate import solve_ivp

from bhwave.conventions import Y_from_Z, Z_from_R


@dataclass
class ScatteringResult:
    omega: float
    x0: float
    x_max: float
    Z: complex | None
    R_boundary: complex | None
    Y: complex
    psi_xmax: complex
    dpsi_xmax: complex
    A_in: complex
    A_out: complex
    R_out: complex
    success: bool
    message: str


def _solve_ivp_with_fallback(*args, method: str, **kwargs):
    try:
        return solve_ivp(*args, method=method, **kwargs)
    except ValueError:
        if method != "RK45":
            return solve_ivp(*args, method="RK45", **kwargs)
        raise


def scattering_reflection(
    omega: float,
    V: Callable[[float], float],
    x0: float,
    x_max: float,
    *,
    Z: complex | None = None,
    R_boundary: complex | None = None,
    psi0: complex = 1.0 + 0j,
    rtol: float = 1e-9,
    atol: float = 1e-11,
    method: str = "DOP853",
    max_step: float | None = None,
) -> ScatteringResult:
    if omega <= 0:
        raise ValueError("omega must be > 0")
    if x_max <= x0:
        raise ValueError("x_max must be greater than x0")
    if (Z is None) == (R_boundary is None):
        raise ValueError("exactly one of Z or R_boundary must be provided")

    if Z is None:
        Z = Z_from_R(R_boundary)
    Y = Y_from_Z(Z, omega)

    def rhs(x: float, y: list[complex]) -> list[complex]:
        psi, dpsi = y
        return [dpsi, (V(x) - omega**2) * psi]

    y0 = [psi0, Y * psi0]
    kwargs = {"rtol": rtol, "atol": atol}
    if max_step is not None:
        kwargs["max_step"] = max_step

    sol = _solve_ivp_with_fallback(rhs, (x0, x_max), y0, method=method, **kwargs)

    psi_xmax = sol.y[0, -1]
    dpsi_xmax = sol.y[1, -1]

    L = x_max - x0
    denom = 2j * omega
    A_in = cmath.exp(1j * omega * L) * (1j * omega * psi_xmax - dpsi_xmax) / denom
    A_out = cmath.exp(-1j * omega * L) * (dpsi_xmax + 1j * omega * psi_xmax) / denom
    R_out = A_out / A_in

    return ScatteringResult(
        omega=omega,
        x0=x0,
        x_max=x_max,
        Z=Z,
        R_boundary=R_boundary,
        Y=Y,
        psi_xmax=psi_xmax,
        dpsi_xmax=dpsi_xmax,
        A_in=A_in,
        A_out=A_out,
        R_out=R_out,
        success=bool(sol.success),
        message=str(sol.message),
    )
