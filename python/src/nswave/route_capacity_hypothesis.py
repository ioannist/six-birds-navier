from __future__ import annotations

import numpy as np


def max_violation_for_p(
    r_steps: np.ndarray,
    p: int,
    *,
    step_offset: int = 0,
) -> float:
    """Return max relative violation for a fixed p and step offset."""
    r_steps = np.asarray(r_steps, dtype=float)
    if r_steps.size == 0:
        return float("nan")
    steps = np.arange(r_steps.size) + step_offset
    eps = 1e-12
    r_model = ((steps + 2) / (steps + 1)) ** p
    ratios = r_steps / np.maximum(r_model, eps)
    return float(np.max(ratios - 1.0))


def fit_smallest_p(
    r_steps: np.ndarray,
    margin_rel: float,
    p_max: int,
    *,
    step_offset: int = 0,
) -> tuple[int | None, float]:
    """Return smallest p meeting r_steps <= r_model*(1+margin), and max violation."""
    r_steps = np.asarray(r_steps, dtype=float)
    if r_steps.size == 0:
        return None, float("nan")
    steps = np.arange(r_steps.size) + step_offset
    eps = 1e-12

    for p in range(p_max + 1):
        r_model = ((steps + 2) / (steps + 1)) ** p
        ratios = r_steps / np.maximum(r_model, eps)
        max_violation = float(np.max(ratios - 1.0))
        if max_violation <= margin_rel:
            return p, max_violation

    return None, max_violation_for_p(r_steps, p_max, step_offset=step_offset)


def _tail_slices(n: int) -> tuple[np.ndarray, np.ndarray]:
    if n <= 0:
        return np.array([], dtype=int), np.array([], dtype=int)
    k_tail = max(4, n // 2)
    idx = np.arange(n)
    return idx, idx[-k_tail:]


def estimate_power_law_exponent(e: np.ndarray, *, eps: float = 1e-12) -> float:
    e = np.asarray(e, dtype=float)
    n = e.size
    if n < 2:
        return float("nan")
    _, tail = _tail_slices(n)
    x = np.log(np.arange(1, n + 1)[tail])
    y = np.log(np.maximum(e[tail], eps))
    A = np.vstack([x, np.ones_like(x)]).T
    coeff, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
    slope = coeff[0]
    return float(-slope)


def tail_ratio_median(e: np.ndarray, *, eps: float = 1e-12) -> float:
    e = np.asarray(e, dtype=float)
    n = e.size
    if n < 2:
        return float("nan")
    _, tail = _tail_slices(n)
    ratios = []
    for i in tail:
        if i + 1 >= n:
            continue
        denom = max(e[i], eps)
        ratios.append(e[i + 1] / denom)
    if not ratios:
        return float("nan")
    return float(np.median(ratios))


def classify_summability(
    e: np.ndarray,
    e_tilde: np.ndarray,
    *,
    min_steps: int = 8,
    eps: float = 1e-12,
) -> dict:
    e = np.asarray(e, dtype=float)
    e_tilde = np.asarray(e_tilde, dtype=float)
    K = e.size
    result = {
        "status": "INCONCLUSIVE",
        "K": int(K),
        "K_tail": int(max(4, K // 2)) if K > 0 else 0,
        "q_est_raw": float("nan"),
        "q_est_tilde": float("nan"),
        "tail_ratio_median_raw": float("nan"),
        "tail_ratio_median_tilde": float("nan"),
    }
    if K < min_steps:
        result["status"] = f"INCONCLUSIVE (need ≥{min_steps} steps; got {K})"
        return result

    result["q_est_raw"] = estimate_power_law_exponent(e, eps=eps)
    result["q_est_tilde"] = estimate_power_law_exponent(e_tilde, eps=eps)
    result["tail_ratio_median_raw"] = tail_ratio_median(e, eps=eps)
    result["tail_ratio_median_tilde"] = tail_ratio_median(e_tilde, eps=eps)

    def _is_decreasing(val: float) -> bool:
        return np.isfinite(val) and val < 0.98

    def _is_plateau(seq: np.ndarray) -> bool:
        if seq.size < 4:
            return False
        increments = np.diff(seq)
        if increments.size < 2:
            return False
        median_inc = float(np.median(increments))
        last_inc = float(increments[-1])
        return median_inc > 0 and last_inc < 0.3 * median_inc

    psum_raw = np.cumsum(e)
    psum_tilde = np.cumsum(e_tilde)
    plateau_raw = _is_plateau(psum_raw)
    plateau_tilde = _is_plateau(psum_tilde)

    pass_raw = (
        result["q_est_raw"] > 1.1
        and _is_decreasing(result["tail_ratio_median_raw"])
        and plateau_raw
    )
    pass_tilde = (
        result["q_est_tilde"] > 1.1
        and _is_decreasing(result["tail_ratio_median_tilde"])
        and plateau_tilde
    )

    fail_raw = result["q_est_raw"] < 0.3 or not _is_decreasing(result["tail_ratio_median_raw"])
    fail_tilde = result["q_est_tilde"] < 0.3 or not _is_decreasing(result["tail_ratio_median_tilde"])

    if pass_raw or pass_tilde:
        result["status"] = "PASS"
    elif fail_raw and fail_tilde:
        result["status"] = "FAIL"
    else:
        result["status"] = "INCONCLUSIVE"
    return result


def assess_mismatch_summability(
    e_route: np.ndarray,
    e_tilde: np.ndarray,
    lambda_direct: np.ndarray,
    *,
    tol_abs: float = 1e-14,
    tol_rel: float = 1e-8,
    min_steps: int = 8,
    eps: float = 1e-12,
) -> dict:
    e_route = np.asarray(e_route, dtype=float)
    e_tilde = np.asarray(e_tilde, dtype=float)
    lambda_direct = np.asarray(lambda_direct, dtype=float)
    scale_candidates = lambda_direct[lambda_direct > 0]
    scale = float(np.median(scale_candidates)) if scale_candidates.size else 0.0
    max_e = float(np.max(e_route)) if e_route.size else 0.0
    thresh = tol_abs + tol_rel * max(scale, 1e-12)

    if max_e <= thresh:
        return {
            "status": "PASS (route mismatch ~0 within tol)",
            "mismatch_trivial": True,
            "K": int(e_route.size),
            "K_tail": int(max(4, e_route.size // 2)) if e_route.size else 0,
            "q_est_raw": None,
            "q_est_tilde": None,
            "tail_ratio_median_raw": None,
            "tail_ratio_median_tilde": None,
            "tol_abs": float(tol_abs),
            "tol_rel": float(tol_rel),
            "scale": float(scale),
        }

    summ = classify_summability(e_route, e_tilde, min_steps=min_steps, eps=eps)
    summ.update(
        {
            "mismatch_trivial": False,
            "tol_abs": float(tol_abs),
            "tol_rel": float(tol_rel),
            "scale": float(scale),
        }
    )
    return summ
