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

from bhwave.barrier import (  # noqa: E402
    barrier_reflection_transmission,
    sample_regge_wheeler_potential,
)


def main() -> None:
    a0 = 0.2
    a2 = 0.3
    b2 = 5.0

    a1_grid = np.linspace(0.0, 1.2, 121)
    b1_grid = np.logspace(np.log10(0.2), np.log10(5.0), 121)

    x, V_rw, x_peak = sample_regge_wheeler_potential(shift_peak_to_zero=True)

    def Vrw_scalar(xq: float) -> float:
        return float(np.interp(xq, x, V_rw))

    x_min = float(x[0])
    x_max = float(x[-1])
    x_ref = 0.0

    x0 = -8.0
    Dx = x_ref - x0

    omega = np.linspace(0.2, 5.0, 81)
    ring_mask = (omega >= 1.0) & (omega <= 3.0)
    low_mask = (omega >= 0.2) & (omega <= 1.0)
    omega_low_target = 0.3
    idx_low = int(np.argmin(np.abs(omega - omega_low_target)))
    omega_low_used = float(omega[idx_low])

    r = np.empty_like(omega, dtype=complex)
    t = np.empty_like(omega, dtype=complex)
    for i, om in enumerate(omega):
        res = barrier_reflection_transmission(
            float(om),
            Vrw_scalar,
            x_min=x_min,
            x_max=x_max,
            x_ref=x_ref,
        )
        r[i] = res.R_barrier
        t[i] = res.T_barrier

    unitarity_resid = np.max(np.abs(np.abs(r) ** 2 + np.abs(t) ** 2 - 1.0))

    a1 = a1_grid[:, None, None]
    b1 = b1_grid[None, :, None]
    om = omega[None, None, :]

    Z = a0 + a1 / (1.0 - 1j * om / b1) + a2 / (1.0 - 1j * om / b2)
    R0 = (1.0 - Z) / (1.0 + Z)

    t_ring = t[ring_mask]
    r_low = r[low_mask]

    D_ring = np.max(np.abs((t_ring**2) * R0[:, :, ring_mask]), axis=2)
    L_proxy = np.abs((a0 + a1_grid + a2) - 1.0)[:, None]
    L_proxy_grid = np.repeat(L_proxy, b1_grid.size, axis=1)
    A_low = 1.0 - np.abs(R0[:, :, idx_low]) ** 2
    G_low = np.max(np.abs(R0[:, :, low_mask] * r_low[None, None, :]), axis=2)

    c1 = D_ring <= 0.05
    c2 = L_proxy_grid <= 0.2
    c3 = A_low >= 0.20
    c4 = G_low <= 0.98

    count = c1.astype(int) + c2.astype(int) + c3.astype(int) + c4.astype(int)

    pct_c1 = 100.0 * float(np.mean(c1))
    pct_c2 = 100.0 * float(np.mean(c2))
    pct_c3 = 100.0 * float(np.mean(c3))
    pct_c4 = 100.0 * float(np.mean(c4))
    pct_all = 100.0 * float(np.mean(count == 4))

    print("BH-18 constraint stacking")
    print(f"  thresholds: D_ring<=0.05, L_proxy<=0.2, A_low>=0.20, G_low<=0.98")
    print(f"  omega_low_used={omega_low_used:.3f}")
    print(f"  max unitarity residual: {unitarity_resid:.3e}")
    print(f"  percent c1 (ring): {pct_c1:.2f}%")
    print(f"  percent c2 (love): {pct_c2:.2f}%")
    print(f"  percent c3 (absorb): {pct_c3:.2f}%")
    print(f"  percent c4 (gain): {pct_c4:.2f}%")
    print(f"  percent all four: {pct_all:.2f}%")

    feasible = np.argwhere(count == 4)
    if feasible.size > 0:
        target = np.array([0.5, 2.0])
        points = np.column_stack((a1_grid[feasible[:, 0]], b1_grid[feasible[:, 1]]))
        d2 = np.sum((points - target) ** 2, axis=1)
        idx_best = feasible[int(np.argmin(d2))]
        i_best, j_best = int(idx_best[0]), int(idx_best[1])
        print(
            "  representative feasible point: "
            f"a1={a1_grid[i_best]:.3f}, b1={b1_grid[j_best]:.3f}"
        )
        print(
            "    metrics: "
            f"D_ring={D_ring[i_best,j_best]:.4f}, "
            f"L_proxy={L_proxy_grid[i_best,j_best]:.4f}, "
            f"A_low={A_low[i_best,j_best]:.4f}, "
            f"G_low={G_low[i_best,j_best]:.4f}"
        )

    B1, A1 = np.meshgrid(b1_grid, a1_grid)

    fig, ax = plt.subplots(figsize=(6.6, 5.2))
    pcm = ax.pcolormesh(B1, A1, count, shading="auto", cmap="viridis", vmin=0, vmax=4)
    fig.colorbar(pcm, ax=ax, label="constraints satisfied")

    ax.contour(B1, A1, D_ring, levels=[0.05], colors="white", linewidths=1.0)
    ax.contour(B1, A1, L_proxy_grid, levels=[0.2], colors="red", linewidths=1.0)
    ax.contour(B1, A1, A_low, levels=[0.20], colors="orange", linewidths=1.0)
    ax.contour(B1, A1, G_low, levels=[0.98], colors="cyan", linewidths=1.0)

    ax.contour(B1, A1, (count == 4).astype(float), levels=[0.5], colors="black", linewidths=2.0)

    ax.set_xscale("log")
    ax.set_xlabel("b1")
    ax.set_ylabel("a1")
    ax.set_title("BH-18 Constraint stacking (toy): Debye Z on RW barrier")

    artifacts_dir = ROOT / "artifacts"
    artifacts_dir.mkdir(parents=True, exist_ok=True)
    plot_path = artifacts_dir / "bh18_constraint_stacking.png"
    fig.tight_layout()
    fig.savefig(plot_path)
    print(f"Saved plot: {plot_path}")


if __name__ == "__main__":
    main()
