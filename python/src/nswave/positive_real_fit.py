from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.optimize import least_squares

Array = np.ndarray


def _softplus(x: Array) -> Array:
    return np.log1p(np.exp(-np.abs(x))) + np.maximum(x, 0)


def _inv_softplus(y: Array) -> Array:
    y = np.asarray(y)
    return np.where(y > 1e-8, np.log(np.expm1(y)), y - 1.0)


def debye_model(omega: Array, a0: float, a: Array, b: Array) -> Array:
    omega = np.asarray(omega)
    H = a0 + 0j
    for ai, bi in zip(a, b):
        H = H + ai / (1.0 - 1j * omega / bi)
    return H


@dataclass
class FitResult:
    a0: float
    a: Array
    b: Array
    success: bool
    rmse: float
    rmse_weighted: float
    min_re_fit: float
    frac_re_nonneg_fit: float


def fit_debye_positive_real(
    omega: Array,
    H: Array,
    *,
    n_modes: int,
    weights: Array | None = None,
    max_nfev: int = 500,
    seed: int = 0,
) -> FitResult:
    omega = np.asarray(omega)
    H = np.asarray(H)
    if weights is None:
        w = np.ones_like(omega, dtype=float)
    else:
        w = np.asarray(weights, dtype=float)
        w = np.sqrt(w / max(np.max(w), 1e-12))

    re_med = float(np.median(np.real(H)))
    a0_init = max(re_med, 0.0)
    a_init = np.full(n_modes, 0.1 * max(1e-6, a0_init + 1.0))
    w_min = float(np.min(omega[omega > 0])) if np.any(omega > 0) else 1.0
    w_max = float(np.max(omega)) if omega.size else 1.0
    b_init = np.logspace(np.log10(w_min + 1e-6), np.log10(w_max + 1e-6), n_modes)

    p0 = np.concatenate(
        [
            np.array([_inv_softplus(a0_init)]),
            _inv_softplus(a_init),
            _inv_softplus(b_init),
        ]
    )

    def unpack(p: Array) -> tuple[float, Array, Array]:
        a0 = float(_softplus(p[0]))
        a = _softplus(p[1: 1 + n_modes])
        b = _softplus(p[1 + n_modes: 1 + 2 * n_modes])
        return a0, a, b

    def residual(p: Array) -> Array:
        a0, a, b = unpack(p)
        H_fit = debye_model(omega, a0, a, b)
        res = (H_fit - H) * w
        return np.concatenate([np.real(res), np.imag(res)])

    rng = np.random.default_rng(seed)
    p0 = p0 + 0.01 * rng.standard_normal(p0.shape)

    result = least_squares(residual, p0, max_nfev=max_nfev)
    a0, a, b = unpack(result.x)
    H_fit = debye_model(omega, a0, a, b)

    res = H_fit - H
    rmse = float(np.sqrt(np.mean(np.abs(res) ** 2)))
    rmse_weighted = float(np.sqrt(np.mean(np.abs(res * w) ** 2)))
    re_vals = np.real(H_fit)
    min_re_fit = float(np.min(re_vals))
    frac_re_nonneg_fit = float(np.mean(re_vals >= -1e-6))

    return FitResult(
        a0=a0,
        a=a,
        b=b,
        success=bool(result.success),
        rmse=rmse,
        rmse_weighted=rmse_weighted,
        min_re_fit=min_re_fit,
        frac_re_nonneg_fit=frac_re_nonneg_fit,
    )
