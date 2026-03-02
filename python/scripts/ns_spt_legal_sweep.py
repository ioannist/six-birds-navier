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

from nswave.deterministic_npz import savez_deterministic  # noqa: E402
from nswave.ns3d_spectral import make_deterministic_u0_hat_3d, simulate_ns3d  # noqa: E402
from nswave.operators import dealias_mask_3d, kgrid_3d, project_div_free_3d  # noqa: E402
from nswave.spt_legal_cert import compute_spt_legal_metrics  # noqa: E402


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


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Stress sweep for NS operational SPT-legal certificate.")
    # Accepted for compatibility with manifest-driven runner interface.
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--resolutions", type=int, nargs="+", default=[32, 48])
    p.add_argument("--nu-list", type=float, nargs="+", default=[0.004, 0.002, 0.001])
    p.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    p.add_argument("--dt", type=float, default=5e-4)
    p.add_argument("--t-max", type=float, default=0.6)
    p.add_argument("--mu", type=float, default=0.0)
    p.add_argument("--alpha", type=float, default=2.0)
    p.add_argument("--k-init", type=int, default=10)
    p.add_argument("--k-force", type=int, default=2)
    p.add_argument("--forcing-amp", type=float, default=1e-2)
    p.add_argument("--record-u-hat-every", type=int, default=30)
    p.add_argument("--window-c", type=float, default=1.0)
    p.add_argument("--energy-threshold", type=float, default=0.95)
    p.add_argument("--kappa-threshold", type=float, default=10.0)
    p.add_argument("--burst-threshold", type=float, default=50.0)
    p.add_argument("--m95-threshold", type=float, default=10.0)
    p.add_argument("--out-dir", default=str(ROOT / "python" / "artifacts"))
    return p.parse_args()


def run_case_metrics(
    *,
    N: int,
    nu: float,
    seed: int,
    dt: float,
    t_max: float,
    mu: float,
    alpha: float,
    k_init: int,
    k_force: int,
    forcing_amp: float,
    record_u_hat_every: int,
    window_c: float,
    energy_threshold: float,
) -> dict:
    L = 2 * np.pi
    u0_hat = make_deterministic_u0_hat_3d(
        N,
        L=L,
        seed=seed,
        k_init=k_init,
        dealias=True,
    )

    if t_max <= 0.0:
        u_hat_hist = np.expand_dims(u0_hat, axis=0)
    else:
        f_hat = make_lowk_divfree_forcing_hat_3d(
            N,
            L=L,
            seed=seed,
            k_force=k_force,
            amp=forcing_amp,
            dealias=True,
        )

        def forcing_hat_fn(
            _t: float, _kx: np.ndarray, _ky: np.ndarray, _kz: np.ndarray
        ) -> np.ndarray:
            return f_hat

        res = simulate_ns3d(
            N=N,
            dt=dt,
            t_max=t_max,
            nu=nu,
            mu=mu,
            alpha=alpha,
            u0_hat=u0_hat,
            forcing_hat_fn=forcing_hat_fn,
            L=L,
            dealias=True,
            record_every=max(1, record_u_hat_every),
            record_u_hat_every=record_u_hat_every,
        )
        if res.u_hat_hist is None:
            raise RuntimeError("u_hat_hist was not recorded")
        u_hat_hist = res.u_hat_hist

    return compute_spt_legal_metrics(
        u_hat_hist,
        N=N,
        L=L,
        window_c=window_c,
        energy_threshold=energy_threshold,
    )


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cfg = {
        "script": Path(__file__).name,
        "resolutions": args.resolutions,
        "nu_list": args.nu_list,
        "seeds": args.seeds,
        "dt": args.dt,
        "t_max": args.t_max,
        "mu": args.mu,
        "alpha": args.alpha,
        "k_init": args.k_init,
        "k_force": args.k_force,
        "forcing_amp": args.forcing_amp,
        "record_u_hat_every": args.record_u_hat_every,
        "window_c": args.window_c,
        "energy_threshold": args.energy_threshold,
        "kappa_threshold": args.kappa_threshold,
        "burst_threshold": args.burst_threshold,
        "m95_threshold": args.m95_threshold,
        "out_dir": str(out_dir),
    }
    cfg_json = json.dumps(cfg, sort_keys=True)
    print(f"FIGURE_CONFIG: {cfg_json}")
    print(f"SPT_LEGAL_SWEEP_CONFIG: {cfg_json}")

    combos = [(N, nu, seed) for N in args.resolutions for nu in args.nu_list for seed in args.seeds]
    total_runs = len(combos)
    print(
        f"runs={total_runs} resolutions={args.resolutions} nu_list={args.nu_list} "
        f"seeds={args.seeds} dt={args.dt} t_max={args.t_max}"
    )

    agg: dict[int, dict[str, list[float]]] = {}
    failures: list[str] = []
    for idx, (N, nu, seed) in enumerate(combos, start=1):
        print(f"[{idx}/{total_runs}] N={N} nu={nu:.4g} seed={seed}")
        try:
            m = run_case_metrics(
                N=N,
                nu=nu,
                seed=seed,
                dt=args.dt,
                t_max=args.t_max,
                mu=args.mu,
                alpha=args.alpha,
                k_init=args.k_init,
                k_force=args.k_force,
                forcing_amp=args.forcing_amp,
                record_u_hat_every=args.record_u_hat_every,
                window_c=args.window_c,
                energy_threshold=args.energy_threshold,
            )
        except Exception as exc:  # pragma: no cover - resilience path
            msg = f"N={N} nu={nu} seed={seed}: {exc}"
            failures.append(msg)
            print(f"  FAIL: {msg}")
            continue

        j_values = m["j_values"]
        kappa = m["kappa_max_over_time"]
        burst = m["burst_ratio"]
        m95 = m["m95"]
        for j, k_val, b_val, m95_val in zip(j_values, kappa, burst, m95):
            jj = int(j)
            if jj not in agg:
                agg[jj] = {"kappa": [], "burst": [], "m95": []}
            agg[jj]["kappa"].append(float(k_val))
            agg[jj]["burst"].append(float(b_val))
            if m95_val is not None:
                agg[jj]["m95"].append(float(m95_val))
        local_max_k = float(np.max(np.asarray(kappa, dtype=float))) if kappa else float("nan")
        local_max_b = float(np.max(np.asarray(burst, dtype=float))) if burst else float("nan")
        local_max_m95 = int(max(m95)) if m95 else 0
        print(
            f"  done: max_kappa={local_max_k:.6e} max_burst={local_max_b:.6e} max_m95={local_max_m95}"
        )

    if not agg:
        print("No successful runs; aborting.")
        return 1

    j_sorted = np.array(sorted(agg.keys()), dtype=int)
    n_cov = np.array([len(agg[int(j)]["kappa"]) for j in j_sorted], dtype=int)
    worst_kappa = np.array([max(agg[int(j)]["kappa"]) for j in j_sorted], dtype=float)
    worst_burst = np.array([max(agg[int(j)]["burst"]) for j in j_sorted], dtype=float)
    worst_m95 = np.array(
        [max(agg[int(j)]["m95"]) if agg[int(j)]["m95"] else np.nan for j in j_sorted],
        dtype=float,
    )

    print("j | runs | worst_kappa_max_over_time | worst_burst_ratio | worst_m95")
    for j, runs, k_val, b_val, m95_val in zip(j_sorted, n_cov, worst_kappa, worst_burst, worst_m95):
        m95_str = "null" if np.isnan(m95_val) else str(int(m95_val))
        print(f"{j:2d} | {runs:4d} | {k_val:25.6e} | {b_val:17.6e} | {m95_str:>9}")

    max_kappa = float(np.max(worst_kappa))
    max_burst = float(np.max(worst_burst))
    plausible = (max_kappa <= args.kappa_threshold) and (max_burst <= args.burst_threshold)
    print(
        "SPT_LEGAL_SWEEP_PLAUSIBLE = "
        f"{plausible} (max_kappa={max_kappa:.6e}, max_burst={max_burst:.6e})"
    )

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    axes[0].plot(j_sorted, worst_kappa, "o-")
    axes[0].axhline(args.kappa_threshold, color="r", linestyle="--", linewidth=1.0)
    axes[0].set_xlabel("dyadic shell j")
    axes[0].set_ylabel("worst kappa_max_over_time")
    axes[0].set_title("SPT-Legal Worst Kappa")
    axes[0].grid(alpha=0.3)

    axes[1].plot(j_sorted, worst_burst, "s-")
    axes[1].axhline(args.burst_threshold, color="r", linestyle="--", linewidth=1.0)
    axes[1].set_xlabel("dyadic shell j")
    axes[1].set_ylabel("worst burst_ratio")
    axes[1].set_title("SPT-Legal Worst Burst Ratio")
    axes[1].grid(alpha=0.3)

    fig.tight_layout()
    out_png = out_dir / "ns_spt_legal_sweep.png"
    plt.savefig(out_png, dpi=150)
    print(f"saved {out_png}")

    out_npz = out_dir / "ns_spt_legal_sweep_data.npz"
    meta_json = np.frombuffer(json.dumps(cfg, sort_keys=True).encode("utf-8"), dtype=np.uint8)
    failures_arr = np.array(failures, dtype=str)
    savez_deterministic(
        out_npz,
        meta_json=meta_json,
        j_sorted=j_sorted,
        n_cov=n_cov,
        worst_kappa=worst_kappa,
        worst_burst=worst_burst,
        worst_m95=worst_m95,
        max_kappa=np.array(max_kappa),
        max_burst=np.array(max_burst),
        plausible=np.array(int(plausible)),
        failures=failures_arr,
    )
    print(f"saved {out_npz}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
