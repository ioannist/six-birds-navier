from __future__ import annotations

import numpy as np


def capacity_L2_ratio(u_hist: np.ndarray, y_hist: np.ndarray, dt: float, eps: float = 1e-12) -> float:
    """
    Estimate the minimal positive-work capacity Λ̂ satisfying:
        ∫ max(Re(<u(t), y(t)>), 0) dt ≤ Λ̂ ∫ ||u(t)||^2 dt
    using discrete time series.
    """
    u = np.asarray(u_hist)
    y = np.asarray(y_hist)
    if u.shape != y.shape:
        raise ValueError("u_hist and y_hist must have the same shape")
    if u.size == 0:
        return 0.0

    u_flat = u.reshape(u.shape[0], -1)
    y_flat = y.reshape(y.shape[0], -1)
    power = np.real(np.sum(np.conj(u_flat) * y_flat, axis=1))
    pos_power = np.maximum(power, 0.0)
    num = np.sum(pos_power) * dt
    denom = np.sum(np.abs(u_flat) ** 2) * dt + eps
    return float(num / denom)


def capacity_sliding_window(
    u_hist: np.ndarray,
    y_hist: np.ndarray,
    dt: float,
    window_len: int,
    hop: int,
    eps: float = 1e-12,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Return (t_mid, Lambda_hat_series) for sliding windows over the time series.
    window_len and hop are in samples.
    """
    if window_len <= 0 or hop <= 0:
        raise ValueError("window_len and hop must be positive")
    u = np.asarray(u_hist)
    y = np.asarray(y_hist)
    if u.shape != y.shape:
        raise ValueError("u_hist and y_hist must have the same shape")
    T = u.shape[0]
    if T < window_len:
        return np.array([], dtype=float), np.array([], dtype=float)

    starts = list(range(0, T - window_len + 1, hop))
    t_mid = np.array([(s + window_len / 2.0) * dt for s in starts], dtype=float)
    lams = []
    for s in starts:
        u_win = u[s : s + window_len]
        y_win = y[s : s + window_len]
        lams.append(capacity_L2_ratio(u_win, y_win, dt, eps=eps))
    return t_mid, np.array(lams, dtype=float)
