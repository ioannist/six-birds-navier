#!/usr/bin/env python3
"""Exact three-mode witness against a shell-only bound for the SGS supplied power.

On the 2*pi torus take k=(4,0,0), l=(-3,1,0), q=k+l=(1,1,0) and
a=(0,1,0), b=(1,3,0), c=(1,-1,0). The datum

    u_eps = a cos(k.x) + b cos(l.x) - eps c sin(q.x)

is smooth, zero mean and divergence free. At j=0 the retained field is
v = G_0 u_eps = -eps c sin(q.x), and the actual supplied SGS output is
y_0 = Delta_0 G_0 (NL(v) - NL(u_eps)) = -(3/2) c sin(q.x). Hence the supplied power
is (3/2) eps |T^3| while the shell energy is eps^2 |T^3|, so their ratio is 3/(2 eps),
whereas the full gradient sup norm stays bounded as eps -> 0.

The script evaluates these quantities with the repository's de-aliased
pseudo-spectral nonlinearity (which returns -NL) and compares them with
the closed forms. Volume-normalized (mean) integrals are reported; sup|grad u| is the maximum over the N^3 grid.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "python" / "src"
sys.path.insert(0, str(SRC))

from nswave.deterministic_npz import savez_deterministic  # noqa: E402
from nswave.dyadic import dyadic_shell_mask  # noqa: E402
from nswave.ns3d_spectral import nonlinear_term_hat_3d  # noqa: E402
from nswave.operators import dealias_mask_3d, kgrid_3d  # noqa: E402

BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
INK = "#0b0b0b"
MUTED = "#52514e"


def witness(N: int, eps_values: np.ndarray) -> dict[str, np.ndarray]:
    x, y, _ = np.meshgrid(*([2 * np.pi * np.arange(N) / N] * 3), indexing="ij")
    a = np.array([0.0, 1.0, 0.0])[:, None, None, None]
    b = np.array([1.0, 3.0, 0.0])[:, None, None, None]
    c = np.array([1.0, -1.0, 0.0])[:, None, None, None]
    high = a * np.cos(4 * x) + b * np.cos(-3 * x + y)
    kx, ky, kz, k2 = kgrid_3d(N)
    shell0 = dyadic_shell_mask(np.sqrt(k2), 0)
    mask = dealias_mask_3d(N)
    kvec = (kx, ky, kz)

    power = np.zeros_like(eps_values)
    shell_energy = np.zeros_like(eps_values)
    grad_sup = np.zeros_like(eps_values)
    output_error = np.zeros_like(eps_values)
    for i, eps in enumerate(eps_values):
        low = -eps * c * np.sin(x + y)
        full = high + low
        full_hat = np.fft.fftn(full, axes=(1, 2, 3))
        low_hat = np.fft.fftn(low, axes=(1, 2, 3))
        # nonlinear_term_hat_3d returns -NL, so this difference is Delta_0 (NL(v) - NL(u)).
        out_hat = (
            nonlinear_term_hat_3d(full_hat, kx, ky, kz, k2, mask=mask)
            - nonlinear_term_hat_3d(low_hat, kx, ky, kz, k2, mask=mask)
        ) * shell0
        out = np.fft.ifftn(out_hat, axes=(1, 2, 3)).real
        output_error[i] = float(np.max(np.abs(out - (-1.5 * c * np.sin(x + y)))))
        power[i] = float(np.mean(np.sum(low * out, axis=0)))
        shell_energy[i] = float(np.mean(np.sum(low * low, axis=0)))
        grads = []
        for comp in range(3):
            for kd in kvec:
                grads.append(np.fft.ifftn(1j * kd * full_hat[comp]).real)
        grad_sup[i] = float(np.max(np.sqrt(np.sum(np.square(grads), axis=0))))
    return {
        "power": power,
        "shell_energy": shell_energy,
        "ratio": power / shell_energy,
        "grad_sup": grad_sup,
        "output_error": output_error,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Exact SGS shell-only bound counterexample")
    parser.add_argument("--N", type=int, default=24)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out-dir", default=str(ROOT / "python" / "artifacts" / "paper"))
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    eps_values = np.logspace(-4, 0, 17)
    cfg = {
        "script": Path(__file__).name,
        "seed": args.seed,
        "N": args.N,
        "dt": 0.0,
        "t_max": 0.0,
        "nu": 0.0,
        "mu": 0.0,
        "alpha": 0.0,
        "shell": 0,
        "eps_min": float(eps_values[0]),
        "eps_max": float(eps_values[-1]),
        "n_eps": int(eps_values.size),
        "out_dir": str(out_dir),
    }
    print(f"FIGURE_CONFIG: {json.dumps(cfg, sort_keys=True)}")

    res = witness(args.N, eps_values)
    max_out_err = float(np.max(res["output_error"]))
    max_ratio_rel_err = float(np.max(np.abs(res["ratio"] * eps_values / 1.5 - 1.0)))
    print("NS-SGS shell-only bound counterexample (j=0, 2*pi torus, mean-normalized integrals)")
    print("eps | supplied power | shell energy | power/shell energy | sup|grad u|")
    for e, p, s, r, g in zip(eps_values, res["power"], res["shell_energy"], res["ratio"], res["grad_sup"]):
        print(f"{e:.3e} | {p:.6e} | {s:.6e} | {r:.6e} | {g:.6f}")
    print(f"max |y_0 - (-(3/2) c sin(q.x))| = {max_out_err:.3e}")
    print(f"max relative error of ratio against 3/(2 eps) = {max_ratio_rel_err:.3e}")
    if max_out_err > 1e-10 or max_ratio_rel_err > 1e-10:
        raise SystemExit("closed form not reproduced")

    plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED})
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.7))
    ax = axes[0]
    ax.loglog(eps_values, res["power"], "o-", color=BLUE, lw=2, ms=5, label=r"supplied power $p_0/|\mathbb{T}^3|$")
    ax.loglog(eps_values, res["shell_energy"], "s-", color=ORANGE, lw=2, ms=5,
              label=r"shell energy $\|u_{\varepsilon,0}\|_2^2/|\mathbb{T}^3|$")
    ax.set_xlabel(r"amplitude $\varepsilon$ of the retained shell")
    ax.set_ylabel("value per unit volume")
    ax.set_title(r"per unit volume: $p_0 = \frac{3}{2}\varepsilon$, $\|u_{\varepsilon,0}\|_2^2=\varepsilon^2$", color=INK, fontsize=10)
    ax.legend(frameon=False)
    ax.grid(alpha=0.25, which="major")

    ax = axes[1]
    ax.loglog(eps_values, res["ratio"] / res["grad_sup"], "o-", color=BLUE, lw=2, ms=5,
              label=r"$p_0 \,/\, (\|\nabla u_\varepsilon\|_\infty \|u_{\varepsilon,0}\|_2^2)$")
    ax.loglog(eps_values, 1.5 / (eps_values * res["grad_sup"]), ":", color=MUTED, lw=1.5,
              label=r"closed form $3/(2\varepsilon\|\nabla u_\varepsilon\|_\infty)$")
    ax.set_xlabel(r"amplitude $\varepsilon$ of the retained shell")
    ax.set_ylabel("ratio a shell-only bound would cap")
    ax.set_title(r"ratio grows like $1/\varepsilon$ while $\|\nabla u_\varepsilon\|_\infty$ stays bounded", color=INK,
                 fontsize=10)
    ax.legend(frameon=False)
    ax.grid(alpha=0.25, which="major")
    gmin, gmax = float(np.min(res["grad_sup"])), float(np.max(res["grad_sup"]))
    ax.text(0.03, 0.05, rf"grid max of $|\nabla u_\varepsilon|$ $\in[{gmin:.2f},\,{gmax:.2f}]$",
            transform=ax.transAxes, color=MUTED, fontsize=9)
    fig.tight_layout()

    png = out_dir / "ns_sgs_shell_counterexample.png"
    fig.savefig(png, dpi=200)
    plt.close(fig)
    print(f"saved {png.resolve()}")

    npz = out_dir / "ns_sgs_shell_counterexample_data.npz"
    meta_json = np.frombuffer(json.dumps(cfg, sort_keys=True).encode("utf-8"), dtype=np.uint8)
    savez_deterministic(
        npz,
        meta_json=meta_json,
        eps=eps_values,
        power=res["power"],
        shell_energy=res["shell_energy"],
        ratio=res["ratio"],
        grad_sup=res["grad_sup"],
        output_error=res["output_error"],
    )
    print(f"saved data: {npz.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
