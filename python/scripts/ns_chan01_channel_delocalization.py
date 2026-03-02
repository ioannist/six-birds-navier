#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "python" / "src"
sys.path.insert(0, str(SRC))

from nswave.channel_diagnostics import (  # noqa: E402
    effective_rank_m95_pooled,
    pca_channels,
)
from nswave.deterministic_npz import savez_deterministic  # noqa: E402
from nswave.dyadic import dyadic_shell_mask  # noqa: E402
from nswave.localization import (  # noqa: E402
    baseline_volume_fraction,
    max_local_energy_fraction,
    window_for_shell,
)
from nswave.ns3d_spectral import make_deterministic_u0_hat_3d, simulate_ns3d  # noqa: E402
from nswave.operators import dealias_mask_3d, kgrid_3d, project_div_free_3d  # noqa: E402


def make_lowk_divfree_forcing_hat_3d(
    N: int,
    *,
    L: float,
    seed: int,
    k_force: int,
    amp: float,
    dealias: bool = True,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    f_phys = rng.standard_normal((3, N, N, N))
    f_hat = np.empty_like(f_phys, dtype=complex)
    for i in range(3):
        f_hat[i] = np.fft.fftn(f_phys[i])
    f_hat[:, 0, 0, 0] = 0.0

    kk = np.fft.fftfreq(N) * N
    kx_i, ky_i, kz_i = np.meshgrid(kk, kk, kk, indexing="ij")
    kmag2 = kx_i * kx_i + ky_i * ky_i + kz_i * kz_i
    f_hat = f_hat * (kmag2 <= k_force**2)

    kx, ky, kz, _ = kgrid_3d(N, L=L)
    f_hat = project_div_free_3d(f_hat, kx, ky, kz)
    if dealias:
        f_hat = f_hat * dealias_mask_3d(N)

    f_phys = [np.fft.ifftn(f_hat[i]).real for i in range(3)]
    rms = np.sqrt(np.mean(f_phys[0] ** 2 + f_phys[1] ** 2 + f_phys[2] ** 2))
    scale = amp / max(rms, 1e-12)
    return f_hat * scale


def shell_real_space(u_hat: np.ndarray, mask: np.ndarray) -> np.ndarray:
    u_shell_hat = u_hat * mask[None, ...]
    out = np.empty_like(u_hat.real)
    for i in range(3):
        out[i] = np.fft.ifftn(u_shell_hat[i]).real
    return out


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="NS-CHAN-01 learned-channel delocalization diagnostic")
    # Accepted for compatibility with manifest-driven runner interface.
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--N", type=int, default=32)
    p.add_argument("--dt", type=float, default=5e-4)
    p.add_argument("--t-max", type=float, default=0.6)
    p.add_argument("--nu", type=float, default=1e-3)
    p.add_argument("--mu", type=float, default=0.0)
    p.add_argument("--alpha", type=float, default=2.0)
    p.add_argument("--k-init", type=int, default=10)
    p.add_argument("--k-force", type=int, default=2)
    p.add_argument("--forcing-amp", type=float, default=1e-2)
    p.add_argument("--record-u-hat-every", type=int, default=20)
    p.add_argument("--burn-in-frac", type=float, default=0.4)
    p.add_argument("--energy-threshold", type=float, default=0.95)
    p.add_argument("--window-c", type=float, default=1.0)
    p.add_argument("--max-k", type=int, default=32)
    p.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    p.add_argument("--out-dir", default=str(ROOT / "python" / "artifacts" / "paper"))
    return p.parse_args()


def main() -> int:
    args = parse_args()
    L = 2 * np.pi
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cfg = {
        "script": Path(__file__).name,
        "N": args.N,
        "dt": args.dt,
        "t_max": args.t_max,
        "nu": args.nu,
        "mu": args.mu,
        "alpha": args.alpha,
        "k_init": args.k_init,
        "k_force": args.k_force,
        "forcing_amp": args.forcing_amp,
        "record_u_hat_every": args.record_u_hat_every,
        "burn_in_frac": args.burn_in_frac,
        "energy_threshold": args.energy_threshold,
        "window_c": args.window_c,
        "max_k": args.max_k,
        "seeds": args.seeds,
        "out_dir": str(out_dir),
    }
    print(f"FIGURE_CONFIG: {json.dumps(cfg, sort_keys=True)}")

    _, _, _, k2 = kgrid_3d(args.N, L=L)
    k_mag = np.sqrt(k2)
    k_max = float(np.max(k_mag))
    j_max = int(np.floor(np.log2(max(k_max, 1.0))))
    j_values = np.arange(1, j_max + 1, dtype=int)
    shell_masks = [dyadic_shell_mask(k_mag, int(j), k0=1.0) for j in j_values]

    hists: list[np.ndarray] = []
    failures: list[str] = []
    for seed in args.seeds:
        try:
            u0_hat = make_deterministic_u0_hat_3d(
                args.N,
                L=L,
                seed=seed,
                k_init=args.k_init,
                dealias=True,
            )
            f_hat = make_lowk_divfree_forcing_hat_3d(
                args.N,
                L=L,
                seed=seed,
                k_force=args.k_force,
                amp=args.forcing_amp,
                dealias=True,
            )

            def forcing_hat_fn(
                _t: float, _kx: np.ndarray, _ky: np.ndarray, _kz: np.ndarray
            ) -> np.ndarray:
                return f_hat

            res = simulate_ns3d(
                N=args.N,
                dt=args.dt,
                t_max=args.t_max,
                nu=args.nu,
                mu=args.mu,
                alpha=args.alpha,
                u0_hat=u0_hat,
                forcing_hat_fn=forcing_hat_fn,
                L=L,
                dealias=True,
                record_every=max(1, args.record_u_hat_every),
                record_u_hat_every=args.record_u_hat_every,
            )
            if res.u_hat_hist is None:
                raise RuntimeError("u_hat_hist was not recorded")
            hist = res.u_hat_hist
            burn_idx = int(np.floor(args.burn_in_frac * hist.shape[0]))
            burn_idx = max(0, min(hist.shape[0] - 1, burn_idx))
            hist = hist[burn_idx:]
            if hist.shape[0] < 1:
                raise RuntimeError("no snapshots after burn-in")
            hists.append(hist)
        except Exception as exc:  # pragma: no cover - resilience path
            failures.append(f"seed={seed}: {exc}")

    if not hists:
        print("NS-CHAN-01 channel delocalization")
        print("no successful runs")
        return 1

    n_feat = 3 * args.N * args.N * args.N
    m95_pooled_vals = np.zeros(len(j_values), dtype=int)
    windows = np.zeros(len(j_values), dtype=int)
    f_vals = np.zeros(len(j_values), dtype=float)
    kappa_max = np.zeros(len(j_values), dtype=float)
    eta_hat = np.zeros(len(j_values), dtype=float)

    for jj, (j, mask) in enumerate(zip(j_values, shell_masks)):
        total_rows = int(sum(h.shape[0] for h in hists))
        X = np.empty((total_rows, n_feat), dtype=np.float32)
        row = 0
        for hist in hists:
            for i in range(hist.shape[0]):
                X[row] = shell_real_space(hist[i], mask).reshape(-1)
                row += 1

        k_init = min(max(1, args.max_k), X.shape[0], X.shape[1])
        m95_pooled = int(effective_rank_m95_pooled(X, energy_threshold=args.energy_threshold))
        svals, V = pca_channels(X, k=k_init)
        _ = svals
        if m95_pooled > V.shape[0]:
            _, V = pca_channels(X, k=m95_pooled)

        window = window_for_shell(args.N, int(j), c=args.window_c)
        base = baseline_volume_fraction(args.N, window)
        kappas: list[float] = []
        for r in range(min(m95_pooled, V.shape[0])):
            psi = V[r].reshape(3, args.N, args.N, args.N)
            conc = max_local_energy_fraction(psi, window)
            kappas.append(float(conc / max(base, 1e-30)))
        kmax = float(max(kappas)) if kappas else 0.0
        eta = float(kmax * m95_pooled * base)

        m95_pooled_vals[jj] = int(m95_pooled)
        windows[jj] = int(window)
        f_vals[jj] = float(base)
        kappa_max[jj] = float(kmax)
        eta_hat[jj] = float(eta)

    print("NS-CHAN-01 channel delocalization")
    print("j | m95_pooled | window | f_j | kappa_max | eta_hat")
    for j, m95, w, f_j, kmax, eta in zip(j_values, m95_pooled_vals, windows, f_vals, kappa_max, eta_hat):
        print(
            f"{int(j):2d} | {int(m95):3d} | {int(w):6d} | {f_j:9.6e} | "
            f"{kmax:9.6e} | {eta:9.6e}"
        )

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    axes[0].plot(j_values, kappa_max, "o-")
    axes[0].set_xlabel("dyadic shell j")
    axes[0].set_ylabel("kappa_max_j")
    axes[0].set_title("Learned Channel Kappa")
    axes[0].grid(alpha=0.3)

    axes[1].plot(j_values, eta_hat, "s-")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("dyadic shell j")
    axes[1].set_ylabel("eta_hat_j")
    axes[1].set_title("Empirical Eta Bound")
    axes[1].grid(alpha=0.3)

    fig.tight_layout()
    png_path = out_dir / "ns_channel_delocalization.png"
    plt.savefig(png_path, dpi=150)
    print(f"saved {png_path}")

    meta_json = np.frombuffer(json.dumps(cfg, sort_keys=True).encode("utf-8"), dtype=np.uint8)
    npz_path = out_dir / "ns_channel_delocalization_data.npz"
    failure_arr = np.array(failures, dtype=str)
    savez_deterministic(
        npz_path,
        meta_json=meta_json,
        j_values=j_values,
        m95_pooled=m95_pooled_vals,
        # Back-compat alias:
        m95=m95_pooled_vals,
        windows=windows,
        baseline_fraction=f_vals,
        kappa_max=kappa_max,
        eta_hat=eta_hat,
        failures=failure_arr,
    )
    print(f"saved {npz_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
