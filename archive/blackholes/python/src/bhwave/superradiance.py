"""Toy superradiance helpers for cavity gain and pole estimates."""

from __future__ import annotations

import cmath
import math

import numpy as np


def round_trip_gain(R0: complex, r: complex) -> float:
    return float(abs(R0 * r))


def pole_estimate_from_G(G: complex, Dx: float, n: int) -> complex:
    if Dx <= 0:
        raise ValueError("Dx must be positive")
    omega_I = (1.0 / (2.0 * Dx)) * math.log(abs(G))
    omega_R = (-cmath.phase(G) + 2.0 * math.pi * n) / (2.0 * Dx)
    return omega_R + 1j * omega_I


def growth_rate_from_gain(gain: float, Dx: float) -> float:
    if Dx <= 0:
        raise ValueError("Dx must be positive")
    if gain <= 0:
        return 0.0
    return max(0.0, float(math.log(gain) / (2.0 * Dx)))
