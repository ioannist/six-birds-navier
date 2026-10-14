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

from bhwave.time_domain import ImpedanceBoundary, simulate_wave_fdtd  # noqa: E402


def main() -> None:
    artifacts_dir = ROOT / "artifacts"
    artifacts_dir.mkdir(parents=True, exist_ok=True)

    V0 = 1.0
    a = 1.0
    x_peak = 0.0
    x_ref = 0.0

    x0 = -8.0
    x_max = 40.0
    Dx = x_ref - x0
    predicted_spacing = 2.0 * Dx
    x_obs = 20.0

    dx = 0.02
    dt = 0.012
    t_max = 120.0
    sponge_start = 30.0
    sponge_gamma_max = 1.0

    boundary = ImpedanceBoundary.from_debye(a0=0.2, a=[0.2], b=[1.0])

    x_src = 10.0
    sigma = 0.8

    def psi_init(x: np.ndarray) -> np.ndarray:
        return np.exp(-((x - x_src) ** 2) / (2.0 * sigma**2))

    def psi_t_init(x: np.ndarray) -> np.ndarray:
        base = psi_init(x)
        return -((x - x_src) / (sigma**2)) * base

    def V_func(x: np.ndarray) -> np.ndarray:
        return V0 / np.cosh(a * (x - x_peak)) ** 2

    print("Config:")
    print(f"  dx={dx}, dt={dt}, CFL={dt/dx:.3f}")
    print(f"  x0={x0}, x_max={x_max}, x_obs={x_obs}")
    print(f"  predicted_spacing={predicted_spacing}")
    print(f"  boundary: a0={boundary.a0}, a={boundary.a.tolist()}, b={boundary.b.tolist()}")
    print("  boundary_update=exact (dt kept at 0.012 for stability)")
    print()

    result = simulate_wave_fdtd(
        V_func=V_func,
        x0=x0,
        x_max=x_max,
        dx=dx,
        dt=dt,
        t_max=t_max,
        boundary=boundary,
        x_obs=x_obs,
        sponge_start=sponge_start,
        sponge_gamma_max=sponge_gamma_max,
        psi_init=psi_init,
        psi_t_init=psi_t_init,
        boundary_update="exact",
    )

    envelope = np.abs(result.psi_obs)
    window = max(1, int(4.0 / dt))
    kernel = np.ones(window) / window
    smooth = np.convolve(envelope, kernel, mode="same")
    threshold = 0.05 * np.max(smooth)
    peaks = []
    min_sep = 0.5 * predicted_spacing
    for i in range(1, len(smooth) - 1):
        if result.t[i] <= 20.0:
            continue
        if smooth[i] < threshold:
            continue
        if smooth[i] >= smooth[i - 1] and smooth[i] >= smooth[i + 1]:
            if not peaks or (result.t[i] - peaks[-1]) >= min_sep:
                peaks.append(result.t[i])

    print("Peak times (first 6, t>20):")
    print("  ", peaks[:6])

    if len(peaks) >= 2:
        spacing = np.median(np.diff(peaks))
        rel_err = abs(spacing - predicted_spacing) / predicted_spacing
        print(f"Median spacing: {spacing:.6g}")
        print(f"Relative error: {rel_err:.6g}")
    else:
        print("Median spacing: n/a")
        print("Relative error: n/a")

    plt.figure(figsize=(7.0, 4.0))
    plt.plot(result.t, envelope, label="|psi_obs|")
    plt.plot(result.t, smooth, label="smoothed", alpha=0.7)
    plt.xlabel("t")
    plt.ylabel("|psi_obs|")
    plt.title("BH-08 time-domain echo demo")
    plt.tight_layout()
    plot_path = artifacts_dir / "bh08_waveform.png"
    plt.savefig(plot_path)
    print(f"Saved plot: {plot_path}")


if __name__ == "__main__":
    main()
