#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "python" / "src"
sys.path.insert(0, str(SRC))

from nswave.ns3d_spectral import make_deterministic_u0_hat_3d, simulate_ns3d  # noqa: E402
from nswave.operators import dealias_mask_3d, kgrid_3d, project_div_free_3d  # noqa: E402
from nswave.spt_legal_cert import (  # noqa: E402
    compute_spt_legal_metrics,
    make_spt_legal_certificate,
    write_cert_json,
)


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
    p = argparse.ArgumentParser(description="Compute operational SPT-legal certificate for NS run.")
    p.add_argument("--N", type=int, default=32)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--nu", type=float, default=1e-3)
    p.add_argument("--mu", type=float, default=0.0)
    p.add_argument("--alpha", type=float, default=2.0)
    p.add_argument("--dt", type=float, default=5e-4)
    p.add_argument("--t-max", type=float, default=0.6)
    p.add_argument("--k-init", type=int, default=10)
    p.add_argument("--k-force", type=int, default=2)
    p.add_argument("--forcing-amp", type=float, default=1e-2)
    p.add_argument("--record-u-hat-every", type=int, default=20)
    p.add_argument("--window-c", type=float, default=1.0)
    p.add_argument("--energy-threshold", type=float, default=0.95)
    p.add_argument("--kappa-threshold", type=float, default=10.0)
    p.add_argument("--burst-threshold", type=float, default=50.0)
    p.add_argument("--require-m95", action="store_true")
    p.add_argument("--m95-threshold", type=float, default=10.0)
    p.add_argument("--no-evolve", action="store_true")
    p.add_argument(
        "--out-json",
        default=str(ROOT / "python" / "artifacts" / "spt_legal_cert.json"),
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    L = 2 * np.pi
    out_json = Path(args.out_json)
    no_evolve = bool(args.no_evolve or args.t_max <= 0.0)

    config = {
        "script": Path(__file__).name,
        "N": args.N,
        "seed": args.seed,
        "nu": args.nu,
        "mu": args.mu,
        "alpha": args.alpha,
        "dt": args.dt,
        "t_max": args.t_max,
        "k_init": args.k_init,
        "k_force": args.k_force,
        "forcing_amp": args.forcing_amp,
        "record_u_hat_every": args.record_u_hat_every,
        "window_c": args.window_c,
        "energy_threshold": args.energy_threshold,
        "kappa_threshold": args.kappa_threshold,
        "burst_threshold": args.burst_threshold,
        "require_m95": bool(args.require_m95),
        "m95_threshold": args.m95_threshold,
        "no_evolve": no_evolve,
        "out_json": str(out_json),
    }
    print(f"SPT_LEGAL_CONFIG: {json.dumps(config, sort_keys=True)}")

    u0_hat = make_deterministic_u0_hat_3d(
        args.N,
        L=L,
        seed=args.seed,
        k_init=args.k_init,
        dealias=True,
    )

    if no_evolve:
        u_hat_hist = np.expand_dims(u0_hat, axis=0)
    else:
        f_hat = make_lowk_divfree_forcing_hat_3d(
            args.N,
            L=L,
            seed=args.seed,
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
        u_hat_hist = res.u_hat_hist

    metrics = compute_spt_legal_metrics(
        u_hat_hist,
        N=args.N,
        L=L,
        window_c=args.window_c,
        energy_threshold=args.energy_threshold,
    )
    thresholds = {
        "kappa_threshold": args.kappa_threshold,
        "burst_threshold": args.burst_threshold,
        "require_m95": bool(args.require_m95),
        "m95_threshold": args.m95_threshold if args.require_m95 else None,
    }
    cert = make_spt_legal_certificate(config=config, metrics=metrics, thresholds=thresholds)
    write_cert_json(cert, out_json)

    max_kappa = cert["summary"]["max_kappa_max_over_time"]
    max_burst = cert["summary"]["max_burst_ratio"]
    print(
        "SPT_LEGAL_STATS: "
        f"max_kappa={max_kappa if max_kappa is not None else 'null'} "
        f"max_burst={max_burst if max_burst is not None else 'null'}"
    )
    print(f"SPT_LEGAL_VERDICT: {cert['verdict']}")
    if cert["reasons"]:
        print("SPT_LEGAL_REASONS:")
        for r in cert["reasons"]:
            print(f"  - {r}")
    print(f"wrote {out_json}")
    return 0 if cert["verdict"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
