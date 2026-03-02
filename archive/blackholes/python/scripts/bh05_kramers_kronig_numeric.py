#!/usr/bin/env python3
from __future__ import annotations

import pathlib
import sys

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

from bhwave.kramers_kronig import kk_reconstruct_re_from_im  # noqa: E402


def debye_sum(omega: np.ndarray, a0: float, a_list: list[float], b_list: list[float]) -> np.ndarray:
    if len(a_list) != len(b_list):
        raise ValueError("a_list and b_list must match")
    Z = np.full_like(omega, a0, dtype=np.complex128)
    for a_k, b_k in zip(a_list, b_list):
        Z = Z + a_k / (1.0 - 1j * omega / b_k)
    return Z


def report_case(name: str, omega: np.ndarray, Z_true: np.ndarray, re_inf: float) -> tuple[float, float]:
    re_true = Z_true.real
    re_pred = kk_reconstruct_re_from_im(omega, Z_true.imag, re_inf=re_inf)
    n = omega.size
    idx0 = int(0.1 * n)
    idx1 = int(0.9 * n)
    err = np.max(np.abs(re_pred[idx0:idx1] - re_true[idx0:idx1]))
    denom = max(1e-12, float(np.max(re_true[idx0:idx1])))
    rel = float(err / denom)
    print(
        f"{name} | N={n} | omega_min={omega[0]:.3e} | omega_max={omega[-1]:.3e} | "
        f"re_inf={re_inf} | max_abs_err_midband={err:.6g} | rel_err_midband={rel:.6g}"
    )
    return err, rel


def main() -> None:
    artifacts_dir = ROOT / "artifacts"
    artifacts_dir.mkdir(parents=True, exist_ok=True)

    omega_grid = np.logspace(-3, 3, 1201)

    cases = [
        ("case1", 0.2, [0.7], [3.0]),
        ("case2", 0.1, [0.4, 0.3, 0.2], [0.8, 3.0, 15.0]),
    ]

    metrics = {}
    for name, a0, a_list, b_list in cases:
        Z_true = debye_sum(omega_grid, a0, a_list, b_list)
        err, rel = report_case(name, omega_grid, Z_true, re_inf=a0)
        metrics[name] = (err, rel)

    # Plot Case 2
    name, a0, a_list, b_list = cases[1]
    Z_true = debye_sum(omega_grid, a0, a_list, b_list)
    re_pred = kk_reconstruct_re_from_im(omega_grid, Z_true.imag, re_inf=a0)
    plt.figure(figsize=(6.0, 4.0))
    plt.plot(omega_grid, Z_true.real, label="Re Z true")
    plt.plot(omega_grid, re_pred, label="Re Z KK", linestyle="--")
    plt.xscale("log")
    plt.xlabel("omega")
    plt.ylabel("Re Z")
    plt.legend()
    plt.tight_layout()
    plot_path = artifacts_dir / "bh05_kk_case2.png"
    plt.savefig(plot_path)
    print(f"Saved plot: {plot_path}")


if __name__ == "__main__":
    main()
