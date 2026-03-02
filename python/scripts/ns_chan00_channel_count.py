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
)
from nswave.deterministic_npz import savez_deterministic  # noqa: E402
from nswave.dyadic import dyadic_shell_mask  # noqa: E402
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
    p = argparse.ArgumentParser(description="NS-CHAN-00 low-rank channel-count diagnostic")
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
    p.add_argument("--window-len", type=int, default=20)
    p.add_argument("--window-hop", type=int, default=5)
    p.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    p.add_argument("--out-dir", default=str(ROOT / "python" / "artifacts" / "paper"))
    return p.parse_args()


def _window_starts(n: int, win_len: int, hop: int) -> list[int]:
    if n <= 0:
        return []
    if n <= win_len:
        return [0]
    starts = list(range(0, n - win_len + 1, hop))
    last = n - win_len
    if starts[-1] != last:
        starts.append(last)
    return starts


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
        "window_len": args.window_len,
        "window_hop": args.window_hop,
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

    # Per-shell aggregate over all seed/window samples.
    m95_samples: dict[int, list[int]] = {int(j): [] for j in j_values}
    m95_pooled_vals: dict[int, int] = {int(j): 0 for j in j_values}
    run_failures: list[str] = []

    hists: list[np.ndarray] = []
    n_feat = 3 * args.N * args.N * args.N
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
            run_failures.append(f"seed={seed}: {exc}")
            continue

        for jj, _mask in enumerate(shell_masks):
            j = int(j_values[jj])
            # Windowed m95 stats remain per-seed for variability diagnostics.
            X = np.empty((hist.shape[0], n_feat), dtype=np.float32)
            for i in range(hist.shape[0]):
                X[i] = shell_real_space(hist[i], shell_masks[jj]).reshape(-1)

            starts = _window_starts(hist.shape[0], args.window_len, args.window_hop)
            for st in starts:
                ed = min(hist.shape[0], st + args.window_len)
                Xw = X[st:ed]
                m95 = effective_rank_m95_pooled(Xw, energy_threshold=args.energy_threshold)
                m95_samples[j].append(int(m95))

    if (not hists) or all(len(v) == 0 for v in m95_samples.values()):
        print("NS-CHAN-00 channel count")
        print("no successful runs")
        return 1

    # Canonical pooled m95 across all post-burn-in snapshots and all seeds.
    for jj, mask in enumerate(shell_masks):
        j = int(j_values[jj])
        total_rows = int(sum(h.shape[0] for h in hists))
        Xp = np.empty((total_rows, n_feat), dtype=np.float32)
        row = 0
        for hist in hists:
            for i in range(hist.shape[0]):
                Xp[row] = shell_real_space(hist[i], mask).reshape(-1)
                row += 1
        m95_pooled_vals[j] = int(
            effective_rank_m95_pooled(Xp, energy_threshold=args.energy_threshold)
        )

    m95_pooled = []
    m95_window_median = []
    m95_window_p25 = []
    m95_window_p75 = []
    m95_window_max = []
    n_windows = []
    n_seeds = []
    for j in j_values:
        vals = np.asarray(m95_samples[int(j)], dtype=float)
        m95_pooled.append(float(m95_pooled_vals[int(j)]))
        if vals.size == 0:
            m95_window_median.append(np.nan)
            m95_window_p25.append(np.nan)
            m95_window_p75.append(np.nan)
            m95_window_max.append(np.nan)
            n_windows.append(0)
        else:
            m95_window_median.append(float(np.median(vals)))
            m95_window_p25.append(float(np.percentile(vals, 25.0)))
            m95_window_p75.append(float(np.percentile(vals, 75.0)))
            m95_window_max.append(float(np.max(vals)))
            n_windows.append(int(vals.size))
        n_seeds.append(int(len(hists)))

    m95_pooled = np.asarray(m95_pooled, dtype=float)
    m95_window_median = np.asarray(m95_window_median, dtype=float)
    m95_window_p25 = np.asarray(m95_window_p25, dtype=float)
    m95_window_p75 = np.asarray(m95_window_p75, dtype=float)
    m95_window_max = np.asarray(m95_window_max, dtype=float)
    n_windows = np.asarray(n_windows, dtype=int)
    n_seeds = np.asarray(n_seeds, dtype=int)

    print("NS-CHAN-00 channel count")
    print(
        "j | m95_pooled | m95_window_median | m95_window_p25 | "
        "m95_window_p75 | m95_window_max | n_windows | n_seeds"
    )
    for j, mp, med, q25, q75, mmax, nw, ns in zip(
        j_values,
        m95_pooled,
        m95_window_median,
        m95_window_p25,
        m95_window_p75,
        m95_window_max,
        n_windows,
        n_seeds,
    ):
        print(
            f"{int(j):2d} | {mp:10.3f} | {med:17.3f} | {q25:14.3f} | {q75:14.3f} | "
            f"{mmax:7.3f} | {int(nw):9d} | {int(ns):7d}"
        )

    fig = plt.figure(figsize=(6.2, 4.4))
    plt.plot(j_values, m95_pooled, "o-", label="m95_pooled (canonical)")
    plt.plot(j_values, m95_window_median, "d-", label="m95_window_median")
    plt.fill_between(j_values, m95_window_p25, m95_window_p75, alpha=0.25, label="window p25-p75")
    plt.plot(j_values, m95_window_max, "s--", linewidth=1.0, label="m95_window_max")
    plt.xlabel("dyadic shell j")
    plt.ylabel("effective channel count m95")
    plt.title("NS-CHAN-00 channel count")
    plt.grid(alpha=0.3)
    plt.legend()
    fig.tight_layout()

    png_path = out_dir / "ns_channel_count.png"
    plt.savefig(png_path, dpi=150)
    print(f"saved {png_path}")

    max_samples = max((len(v) for v in m95_samples.values()), default=0)
    m95_sample_matrix = np.full((len(j_values), max_samples), -1, dtype=int)
    for idx, j in enumerate(j_values):
        vals = np.asarray(m95_samples[int(j)], dtype=int)
        if vals.size:
            m95_sample_matrix[idx, : vals.size] = vals

    meta_json = np.frombuffer(json.dumps(cfg, sort_keys=True).encode("utf-8"), dtype=np.uint8)
    npz_path = out_dir / "ns_channel_count_data.npz"
    failure_arr = np.array(run_failures, dtype=str)
    savez_deterministic(
        npz_path,
        meta_json=meta_json,
        j_values=j_values,
        m95_pooled=m95_pooled,
        m95_window_median=m95_window_median,
        m95_window_p25=m95_window_p25,
        m95_window_p75=m95_window_p75,
        m95_window_max=m95_window_max,
        # Back-compat aliases:
        m95_median=m95_window_median,
        m95_p25=m95_window_p25,
        m95_p75=m95_window_p75,
        m95_max=m95_window_max,
        n_windows=n_windows,
        n_seeds=n_seeds,
        m95_samples=m95_sample_matrix,
        failures=failure_arr,
    )
    print(f"saved {npz_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
