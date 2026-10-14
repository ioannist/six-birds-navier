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

from bhwave.kramers_kronig import kk_reconstruct_re_from_im, kk_static_re0_from_im  # noqa: E402


def Z_debye_sum(omega: np.ndarray, a0: float, a_list: list[float], b_list: list[float]) -> np.ndarray:
    if len(a_list) != len(b_list):
        raise ValueError("a_list and b_list must match")
    Z = np.full_like(omega, a0, dtype=np.complex128)
    for a_k, b_k in zip(a_list, b_list):
        Z = Z + a_k / (1.0 - 1j * omega / b_k)
    return Z


def main() -> None:
    artifacts_dir = ROOT / "artifacts"
    artifacts_dir.mkdir(parents=True, exist_ok=True)

    omega = np.logspace(-4, 4, 4001)

    cases = [
        ("case1", 0.2, [0.7], [3.0]),
        ("case2", 0.1, [0.4, 0.3, 0.2], [0.8, 3.0, 15.0]),
        ("case3", 0.1, [0.2, 0.15, 0.1], [0.8, 3.0, 15.0]),
    ]

    print("case | ReZ0_true | ReZ0_pred | abs_err | delta_true | delta_pred")
    for name, a0, a_list, b_list in cases:
        Z = Z_debye_sum(omega, a0, a_list, b_list)
        re_inf = a0
        ReZ0_true = a0 + sum(a_list)
        ReZ0_pred = kk_static_re0_from_im(omega, Z.imag, re_inf=re_inf)
        abs_err = abs(ReZ0_pred - ReZ0_true)
        delta_true = sum(a_list)
        delta_pred = ReZ0_pred - a0
        print(
            f"{name} | {ReZ0_true:.6g} | {ReZ0_pred:.6g} | {abs_err:.6g} | "
            f"{delta_true:.6g} | {delta_pred:.6g}"
        )

    print(
        "Static conservative proxy (ReZ0) is reconstructed from dissipative spectrum ImZ/omega via KK."
    )

    # Plot for case2
    _, a0, a_list, b_list = cases[1]
    Z = Z_debye_sum(omega, a0, a_list, b_list)
    re_pred_grid = kk_reconstruct_re_from_im(omega, Z.imag, re_inf=a0)

    plt.figure(figsize=(6.0, 4.0))
    plt.plot(omega, Z.real, label="Re Z true")
    plt.plot(omega, re_pred_grid, label="Re Z KK", linestyle="--")
    plt.xscale("log")
    plt.xlabel("omega")
    plt.ylabel("Re Z")
    plt.legend()
    plt.tight_layout()
    plot_path = artifacts_dir / "bh10_cross_channel.png"
    plt.savefig(plot_path)
    print(f"Saved plot: {plot_path}")


if __name__ == "__main__":
    main()
