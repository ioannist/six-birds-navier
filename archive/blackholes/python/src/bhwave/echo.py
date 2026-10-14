"""Echo/cavity transfer helpers for barrier + inner boundary models."""

from __future__ import annotations

import cmath
import math
import warnings


def _check_inputs(omega: float, Dx: float) -> None:
    if omega <= 0:
        raise ValueError("omega must be > 0")
    if not math.isfinite(Dx):
        raise ValueError("Dx must be finite")


def cavity_H(omega: float, R0: complex, r: complex, Dx: float) -> complex:
    _check_inputs(omega, Dx)
    denom = 1.0 - R0 * r * cmath.exp(2j * omega * Dx)
    if abs(denom) < 1e-12:
        warnings.warn("cavity_H denominator near zero", RuntimeWarning)
    return 1.0 / denom


def predicted_R_out_x0(
    omega: float,
    R0: complex,
    r: complex,
    t: complex,
    Dx: float,
) -> complex:
    _check_inputs(omega, Dx)
    denom = 1.0 - R0 * r * cmath.exp(2j * omega * Dx)
    if abs(denom) < 1e-12:
        warnings.warn("predicted_R_out_x0 denominator near zero", RuntimeWarning)
    return r * cmath.exp(-2j * omega * Dx) + (t * t * R0) / denom
