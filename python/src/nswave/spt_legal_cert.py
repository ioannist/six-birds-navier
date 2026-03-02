"""SPT-legal anti-localization certificate helpers."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from nswave.dyadic import dyadic_shell_mask
from nswave.localization import (
    localization_factors,
    max_local_energy_fraction,
    window_for_shell,
)
from nswave.operators import kgrid_3d

Array = np.ndarray


def _shell_real_space(u_hat: Array, mask: Array) -> Array:
    u_shell_hat = u_hat * mask[None, ...]
    out = np.empty_like(u_hat.real)
    for i in range(3):
        out[i] = np.fft.ifftn(u_shell_hat[i]).real
    return out


def _m95_from_snapshot_matrix(X: Array, energy_threshold: float) -> int:
    x_norm = float(np.linalg.norm(X))
    if x_norm <= 1e-14:
        return 0
    svals = np.linalg.svd(X, full_matrices=False, compute_uv=False)
    s2 = svals * svals
    total = float(np.sum(s2))
    if total <= 0.0:
        return 0
    cum = np.cumsum(s2) / total
    return int(np.searchsorted(cum, energy_threshold) + 1)


def _as_float_list(x: Array) -> list[float]:
    return [float(v) for v in np.asarray(x, dtype=float).tolist()]


def _as_int_list(x: Array) -> list[int]:
    return [int(v) for v in np.asarray(x, dtype=int).tolist()]


def compute_spt_legal_metrics(
    u_hat_hist: Array,
    *,
    N: int,
    L: float,
    window_c: float,
    energy_threshold: float,
) -> dict[str, Any]:
    """
    Compute per-shell anti-localization metrics from u_hat snapshots.

    Returns a pure-Python dict suitable for deterministic JSON serialization.
    """
    hist = np.asarray(u_hat_hist)
    if hist.ndim != 5 or hist.shape[1] != 3:
        raise ValueError("u_hat_hist must have shape (T, 3, N, N, N)")
    if hist.shape[2] != N or hist.shape[3] != N or hist.shape[4] != N:
        raise ValueError("u_hat_hist shape does not match N")
    if window_c <= 0:
        raise ValueError("window_c must be positive")
    if not (0.0 < energy_threshold <= 1.0):
        raise ValueError("energy_threshold must lie in (0, 1]")

    n_snap = int(hist.shape[0])
    _, _, _, k2 = kgrid_3d(N, L=L)
    k_mag = np.sqrt(k2)
    k_max = float(np.max(k_mag))
    j_max = int(np.floor(np.log2(max(k_max, 1.0))))
    j_values = np.arange(1, j_max + 1, dtype=int)

    if j_values.size == 0:
        return {
            "N": int(N),
            "L": float(L),
            "window_c": float(window_c),
            "energy_threshold": float(energy_threshold),
            "n_snapshots": n_snap,
            "j_values": [],
            "windows": [],
            "baseline_fraction": [],
            "mean_conc": [],
            "max_conc": [],
            "L_mean": [],
            "L_max": [],
            "kappa_mean": [],
            "kappa_max_over_time": [],
            "burst_ratio": [],
            "m95": [],
        }

    shell_masks = [dyadic_shell_mask(k_mag, int(j), k0=1.0) for j in j_values]
    windows = [window_for_shell(N, int(j), c=window_c) for j in j_values]
    n_shell = len(j_values)
    n_feat = 3 * N * N * N

    conc = np.zeros((n_snap, n_shell), dtype=float)
    m95 = np.zeros(n_shell, dtype=int)

    for jj, (mask, window) in enumerate(zip(shell_masks, windows)):
        X = np.empty((n_snap, n_feat), dtype=np.float32)
        for i in range(n_snap):
            u_shell = _shell_real_space(hist[i], mask)
            conc[i, jj] = max_local_energy_fraction(u_shell, window)
            X[i] = u_shell.reshape(-1)
        m95[jj] = _m95_from_snapshot_matrix(X, energy_threshold=energy_threshold)

    mean_conc = np.mean(conc, axis=0)
    max_conc = np.max(conc, axis=0)
    f_j, L_mean, kappa_mean = localization_factors(
        mean_conc,
        N=N,
        windows=np.array(windows, dtype=int),
        j_values=j_values,
    )
    _, L_max, kappa_max = localization_factors(
        max_conc,
        N=N,
        windows=np.array(windows, dtype=int),
        j_values=j_values,
    )
    burst_ratio = max_conc / np.maximum(mean_conc, 1e-30)

    return {
        "N": int(N),
        "L": float(L),
        "window_c": float(window_c),
        "energy_threshold": float(energy_threshold),
        "n_snapshots": n_snap,
        "j_values": _as_int_list(j_values),
        "windows": [int(w) for w in windows],
        "baseline_fraction": _as_float_list(f_j),
        "mean_conc": _as_float_list(mean_conc),
        "max_conc": _as_float_list(max_conc),
        "L_mean": _as_float_list(L_mean),
        "L_max": _as_float_list(L_max),
        "kappa_mean": _as_float_list(kappa_mean),
        "kappa_max_over_time": _as_float_list(kappa_max),
        "burst_ratio": _as_float_list(burst_ratio),
        "m95": _as_int_list(m95),
    }


def certificate_verdict(cert: dict[str, Any]) -> tuple[bool, list[str]]:
    thresholds = cert.get("thresholds", {})
    metrics = cert.get("metrics", {})

    kappa_vals = np.asarray(metrics.get("kappa_max_over_time", []), dtype=float)
    burst_vals = np.asarray(metrics.get("burst_ratio", []), dtype=float)
    m95_vals = [int(v) for v in metrics.get("m95", []) if v is not None]

    reasons: list[str] = []
    if kappa_vals.size == 0:
        reasons.append("no shell metrics available (empty kappa_max_over_time)")
    if burst_vals.size == 0:
        reasons.append("no shell metrics available (empty burst_ratio)")

    kappa_threshold = float(thresholds.get("kappa_threshold", 10.0))
    burst_threshold = float(thresholds.get("burst_threshold", 50.0))
    require_m95 = bool(thresholds.get("require_m95", False))
    m95_threshold_raw = thresholds.get("m95_threshold")
    m95_threshold = None if m95_threshold_raw is None else float(m95_threshold_raw)

    if kappa_vals.size:
        max_kappa = float(np.max(kappa_vals))
        if max_kappa > kappa_threshold:
            reasons.append(
                f"max_kappa_max_over_time={max_kappa:.6e} exceeds kappa_threshold={kappa_threshold:.6e}"
            )
    if burst_vals.size:
        max_burst = float(np.max(burst_vals))
        if max_burst > burst_threshold:
            reasons.append(
                f"max_burst_ratio={max_burst:.6e} exceeds burst_threshold={burst_threshold:.6e}"
            )

    if require_m95:
        if not m95_vals:
            reasons.append("require_m95=True but m95 values are unavailable")
        elif m95_threshold is None:
            reasons.append("require_m95=True but m95_threshold is not set")
        else:
            max_m95 = int(max(m95_vals))
            if max_m95 > m95_threshold:
                reasons.append(
                    f"max_m95={max_m95} exceeds m95_threshold={int(m95_threshold)}"
                )

    return (len(reasons) == 0), reasons


def make_spt_legal_certificate(
    config: dict[str, Any],
    metrics: dict[str, Any],
    thresholds: dict[str, Any],
) -> dict[str, Any]:
    kappa_vals = np.asarray(metrics.get("kappa_max_over_time", []), dtype=float)
    burst_vals = np.asarray(metrics.get("burst_ratio", []), dtype=float)
    m95_vals = [int(v) for v in metrics.get("m95", []) if v is not None]

    summary = {
        "n_shells": int(len(metrics.get("j_values", []))),
        "n_snapshots": int(metrics.get("n_snapshots", 0)),
        "max_kappa_max_over_time": float(np.max(kappa_vals)) if kappa_vals.size else None,
        "max_burst_ratio": float(np.max(burst_vals)) if burst_vals.size else None,
        "max_m95": int(max(m95_vals)) if m95_vals else None,
    }

    cert = {
        "schema_version": "ns_spt_legal_cert.v1",
        "config": config,
        "thresholds": thresholds,
        "summary": summary,
        "metrics": metrics,
    }
    ok, reasons = certificate_verdict(cert)
    cert["verdict"] = "PASS" if ok else "FAIL"
    cert["reasons"] = reasons
    return cert


def write_cert_json(cert: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(cert, sort_keys=True, indent=2)
    path.write_text(text + "\n", encoding="utf-8")
