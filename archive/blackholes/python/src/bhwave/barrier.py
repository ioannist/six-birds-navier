"""Barrier-only scattering for toy black-hole potentials."""

from __future__ import annotations

from dataclasses import dataclass
import cmath
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import lambertw


def poeschl_teller_V(x: float, V0: float, a: float, x_peak: float = 0.0) -> float:
    return float(V0 / np.cosh(a * (x - x_peak)) ** 2)


def rstar_of_r(r: float | np.ndarray, M: float = 1.0) -> float | np.ndarray:
    r_arr = np.asarray(r, dtype=float)
    if np.any(r_arr <= 2.0 * M):
        raise ValueError("r must be greater than 2M")
    rstar = r_arr + 2.0 * M * np.log(r_arr / (2.0 * M) - 1.0)
    if np.ndim(r_arr) == 0:
        return float(rstar)
    return rstar


def r_of_rstar(rstar: float | np.ndarray, M: float = 1.0) -> float | np.ndarray:
    rstar_arr = np.asarray(rstar, dtype=float)
    arg = np.exp(rstar_arr / (2.0 * M) - 1.0)
    u = lambertw(arg)
    r = 2.0 * M * (1.0 + u)
    if np.ndim(rstar_arr) == 0:
        return float(np.real(r))
    return np.real(r)


def regge_wheeler_V_r(r: float | np.ndarray, M: float = 1.0, ell: int = 2) -> float | np.ndarray:
    r_arr = np.asarray(r, dtype=float)
    f = 1.0 - 2.0 * M / r_arr
    V = f * (ell * (ell + 1) / r_arr**2 - 6.0 * M / r_arr**3)
    if np.ndim(r_arr) == 0:
        return float(V)
    return V


def regge_wheeler_V_rstar(x: float | np.ndarray, M: float = 1.0, ell: int = 2) -> float | np.ndarray:
    r = r_of_rstar(x, M)
    return regge_wheeler_V_r(r, M, ell)


def fit_poeschl_teller_from_samples(x: np.ndarray, V: np.ndarray) -> tuple[float, float, float]:
    x_arr = np.asarray(x, dtype=float)
    V_arr = np.asarray(V, dtype=float)
    if x_arr.ndim != 1 or V_arr.ndim != 1 or x_arr.shape != V_arr.shape:
        raise ValueError("x and V must be 1D arrays of the same shape")
    if x_arr.size < 3:
        raise ValueError("need at least 3 samples to estimate curvature")
    i_peak = int(np.argmax(V_arr))
    if i_peak == 0 or i_peak == x_arr.size - 1:
        raise ValueError("peak index must have neighbors")
    V0 = float(V_arr[i_peak])
    if V0 <= 0:
        raise ValueError("peak height must be positive")
    dx = float(x_arr[i_peak + 1] - x_arr[i_peak])
    if dx == 0.0:
        raise ValueError("x grid spacing must be nonzero")
    Vpp = (V_arr[i_peak + 1] - 2.0 * V_arr[i_peak] + V_arr[i_peak - 1]) / (dx * dx)
    if Vpp >= 0:
        raise ValueError("peak curvature must be negative")
    a = float(np.sqrt(-Vpp / (2.0 * V0)))
    return float(x_arr[i_peak]), V0, a


def sample_regge_wheeler_potential(
    *,
    M: float = 1.0,
    ell: int = 2,
    x_min: float = -80.0,
    x_max: float = 140.0,
    n: int = 22001,
    shift_peak_to_zero: bool = True,
) -> tuple[np.ndarray, np.ndarray, float]:
    """
    Return (x, V(x), x_peak_unshifted) on a uniform tortoise grid.
    If shift_peak_to_zero, x is shifted so argmax(V) is at x=0.
    """
    if n < 3:
        raise ValueError("n must be at least 3")
    if x_max <= x_min:
        raise ValueError("x_max must be greater than x_min")
    x_un = np.linspace(x_min, x_max, n)
    V_un = regge_wheeler_V_rstar(x_un, M=M, ell=ell)
    i_peak = int(np.argmax(V_un))
    x_peak_unshifted = float(x_un[i_peak])
    if shift_peak_to_zero:
        x = x_un - x_peak_unshifted
    else:
        x = x_un
    return x, np.asarray(V_un, dtype=float), x_peak_unshifted


@dataclass
class BarrierScatteringResult:
    omega: float
    x_min: float
    x_max: float
    x_ref: float
    psi_xmin: complex
    dpsi_xmin: complex
    A_in: complex
    A_out: complex
    R_barrier: complex
    T_barrier: complex
    success: bool
    message: str


def _solve_ivp_with_fallback(*args, method: str, **kwargs):
    try:
        return solve_ivp(*args, method=method, **kwargs)
    except ValueError:
        if method != "RK45":
            return solve_ivp(*args, method="RK45", **kwargs)
        raise


def barrier_reflection_transmission(
    omega: float,
    V: Callable[[float], float],
    *,
    x_min: float,
    x_max: float,
    x_ref: float = 0.0,
    rtol: float = 1e-9,
    atol: float = 1e-11,
    method: str = "DOP853",
    max_step: float | None = None,
) -> BarrierScatteringResult:
    if omega <= 0:
        raise ValueError("omega must be > 0")
    if x_max <= x_min:
        raise ValueError("x_max must be greater than x_min")

    yR = x_max - x_ref
    psi0 = cmath.exp(1j * omega * yR)
    dpsi0 = 1j * omega * psi0

    def rhs(x: float, y: list[complex]) -> list[complex]:
        psi, dpsi = y
        return [dpsi, (V(x) - omega**2) * psi]

    y0 = [psi0, dpsi0]
    kwargs = {"rtol": rtol, "atol": atol}
    if max_step is not None:
        kwargs["max_step"] = max_step

    sol = _solve_ivp_with_fallback(rhs, (x_max, x_min), y0, method=method, **kwargs)

    psi_xmin = sol.y[0, -1]
    dpsi_xmin = sol.y[1, -1]

    yL = x_min - x_ref
    Eplus = cmath.exp(1j * omega * yL)
    Eminus = cmath.exp(-1j * omega * yL)
    denom = 1j * omega
    A_in = (psi_xmin + dpsi_xmin / denom) / (2.0 * Eplus)
    A_out = (psi_xmin - dpsi_xmin / denom) / (2.0 * Eminus)

    R_barrier = A_out / A_in
    T_barrier = 1.0 / A_in

    return BarrierScatteringResult(
        omega=omega,
        x_min=x_min,
        x_max=x_max,
        x_ref=x_ref,
        psi_xmin=psi_xmin,
        dpsi_xmin=dpsi_xmin,
        A_in=A_in,
        A_out=A_out,
        R_barrier=R_barrier,
        T_barrier=T_barrier,
        success=bool(sol.success),
        message=str(sol.message),
    )
