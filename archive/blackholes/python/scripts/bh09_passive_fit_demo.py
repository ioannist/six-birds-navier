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

from bhwave.barrier import barrier_reflection_transmission, poeschl_teller_V  # noqa: E402
from bhwave.conventions import R_from_Z  # noqa: E402
from bhwave.inference import (
    R_pred_grid,
    Z_debye_k1,
    estimate_Rpred_from_waveform,
    fit_debye_k1_from_Rpred,
    gaussian_source,
    passivity_residuals,
    synthesize_waveform_irfft,
)  # noqa: E402


def main() -> None:
    artifacts_dir = ROOT / "artifacts"
    artifacts_dir.mkdir(parents=True, exist_ok=True)

    N = 512
    dt = 0.25
    omega_pos = 2.0 * np.pi * np.fft.rfftfreq(N, dt)

    V0 = 1.0
    a = 1.0
    x_peak = 0.0
    x_ref = 0.0
    x_min = -12.0
    x_max = 12.0
    x0 = -8.0
    Dx = x_ref - x0

    def V(x: float) -> float:
        return poeschl_teller_V(x, V0, a, x_peak)

    r = np.zeros_like(omega_pos, dtype=complex)
    t = np.zeros_like(omega_pos, dtype=complex)
    for i, omega in enumerate(omega_pos):
        if omega <= 0:
            r[i] = 1.0 + 0j
            t[i] = 0.0 + 0j
            continue
        res = barrier_reflection_transmission(
            float(omega),
            V,
            x_min=x_min,
            x_max=x_max,
            x_ref=x_ref,
        )
        r[i] = res.R_barrier
        t[i] = res.T_barrier

    a0_true = 0.2
    a1_true = 0.5
    b1_true = 2.0

    Z_true = Z_debye_k1(omega_pos, a0_true, a1_true, b1_true)
    R0_true = np.vectorize(R_from_Z)(Z_true)
    Rpred_true = R_pred_grid(omega_pos, R0_true, r, t, Dx)

    S_pos = gaussian_source(omega_pos, omega0=1.5, sigma=0.5)
    S_pos[0] = 0.0

    tgrid, y, _ = synthesize_waveform_irfft(N=N, dt=dt, S_pos=S_pos, Rpred_pos=Rpred_true)

    Rpred_est = estimate_Rpred_from_waveform(waveform=y, S_pos=S_pos)

    p0 = (0.15, 0.3, 1.0)
    weights = np.abs(S_pos)
    if np.max(weights) > 0:
        weights = weights / np.max(weights)

    fit = fit_debye_k1_from_Rpred(
        omega_pos=omega_pos,
        Rpred_target=Rpred_est,
        r=r,
        t=t,
        Dx=Dx,
        p0=p0,
        weight=weights,
    )

    a0_fit = fit["a0"]
    a1_fit = fit["a1"]
    b1_fit = fit["b1"]

    Z_fit = Z_debye_k1(omega_pos, a0_fit, a1_fit, b1_fit)
    R0_fit = np.vectorize(R_from_Z)(Z_fit)
    Rpred_fit = R_pred_grid(omega_pos, R0_fit, r, t, Dx)

    resid_R, resid_Z = passivity_residuals(omega_pos, Z_fit, R0_fit)

    mask = np.abs(S_pos) > 1e-3 * np.max(np.abs(S_pos))
    max_abs_err = np.max(np.abs(Rpred_fit[mask] - Rpred_est[mask]))

    print("Config:")
    print(f"  N={N}, dt={dt}, Dx={Dx}")
    print(f"  barrier: V0={V0}, a={a}, x_peak={x_peak}, x_min={x_min}, x_max={x_max}")
    print(f"  true params: a0={a0_true}, a1={a1_true}, b1={b1_true}")
    print(f"  init guess: a0={p0[0]}, a1={p0[1]}, b1={p0[2]}")
    print()

    def rel_err(est: float, true: float) -> float:
        return abs(est - true) / max(1e-12, abs(true))

    print("Fit results:")
    print(f"  recovered: a0={a0_fit:.6g}, a1={a1_fit:.6g}, b1={b1_fit:.6g}")
    print(
        "  rel errors: a0={:.6g}, a1={:.6g}, b1={:.6g}".format(
            rel_err(a0_fit, a0_true),
            rel_err(a1_fit, a1_true),
            rel_err(b1_fit, b1_true),
        )
    )
    print(f"  passivity residuals: resid_R={resid_R:.6g}, resid_Z={resid_Z:.6g}")
    print(f"  max_abs_err_Rpred(masked)={max_abs_err:.6g}")

    _, y_fit, _ = synthesize_waveform_irfft(N=N, dt=dt, S_pos=S_pos, Rpred_pos=Rpred_fit)

    plt.figure(figsize=(7.0, 4.0))
    plt.plot(tgrid, y, label="waveform")
    plt.plot(tgrid, y_fit, label="fit", linestyle="--")
    plt.xlabel("t")
    plt.ylabel("psi")
    plt.title("BH-09 passive fit demo")
    plt.legend()
    plt.tight_layout()
    plot_path = artifacts_dir / "bh09_fit_waveform.png"
    plt.savefig(plot_path)
    print(f"Saved plot: {plot_path}")


if __name__ == "__main__":
    main()
