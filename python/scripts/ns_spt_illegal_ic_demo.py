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
from nswave.illegal_initial_conditions import (  # noqa: E402
    bandlimit_to_shell,
    make_localized_blob_u0_hat,
    rescale_energy,
)
from nswave.ns3d_spectral import make_deterministic_u0_hat_3d  # noqa: E402
from nswave.spt_legal_cert import (  # noqa: E402
    make_spt_legal_certificate,
    compute_spt_legal_metrics,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Demo: intentionally SPT-illegal no-evolve initial condition.")
    p.add_argument("--N", type=int, default=32)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--k-init", type=int, default=10)
    p.add_argument("--j-target", type=int, default=3)
    p.add_argument("--sigma", type=float, default=None)
    p.add_argument("--center-x", type=float, default=0.23)
    p.add_argument("--center-y", type=float, default=0.41)
    p.add_argument("--center-z", type=float, default=0.67)
    p.add_argument("--illegal-dealias", action="store_true")
    p.add_argument("--window-c", type=float, default=1.0)
    p.add_argument("--energy-threshold", type=float, default=0.95)
    p.add_argument("--kappa-threshold", type=float, default=10.0)
    p.add_argument("--burst-threshold", type=float, default=50.0)
    p.add_argument(
        "--out-dir",
        default=str(ROOT / "python" / "artifacts" / "paper"),
    )
    p.add_argument("--no-evolve", action="store_true")
    return p.parse_args()


def _cert_from_u0(
    u0_hat: np.ndarray,
    *,
    N: int,
    window_c: float,
    energy_threshold: float,
    kappa_threshold: float,
    burst_threshold: float,
    config: dict,
) -> dict:
    u_hist = np.expand_dims(u0_hat, axis=0)
    metrics = compute_spt_legal_metrics(
        u_hist,
        N=N,
        L=2 * np.pi,
        window_c=window_c,
        energy_threshold=energy_threshold,
    )
    thresholds = {
        "kappa_threshold": float(kappa_threshold),
        "burst_threshold": float(burst_threshold),
        "require_m95": False,
        "m95_threshold": None,
    }
    return make_spt_legal_certificate(config=config, metrics=metrics, thresholds=thresholds)


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sigma = float(args.sigma) if args.sigma is not None else (0.8 / float(args.N))
    center = (float(args.center_x), float(args.center_y), float(args.center_z))
    no_evolve = True

    cfg = {
        "script": Path(__file__).name,
        "N": args.N,
        "seed": args.seed,
        "k_init": args.k_init,
        "j_target": args.j_target,
        "sigma": sigma,
        "center_xyz": [center[0], center[1], center[2]],
        "illegal_dealias": bool(args.illegal_dealias),
        "window_c": args.window_c,
        "energy_threshold": args.energy_threshold,
        "kappa_threshold": args.kappa_threshold,
        "burst_threshold": args.burst_threshold,
        "no_evolve": no_evolve,
        "out_dir": str(out_dir),
    }
    print(f"FIGURE_CONFIG: {json.dumps(cfg, sort_keys=True)}")

    legal_u0_hat = make_deterministic_u0_hat_3d(
        args.N,
        L=2 * np.pi,
        seed=args.seed,
        k_init=args.k_init,
        dealias=True,
    )

    illegal_raw = make_localized_blob_u0_hat(
        args.N,
        center_xyz=center,
        sigma=sigma,
        seed=args.seed,
    )
    illegal_shell = bandlimit_to_shell(
        illegal_raw,
        j=args.j_target,
        dealias=bool(args.illegal_dealias),
    )
    # Normalize illegal energy against legal baseline for apples-to-apples thresholding.
    legal_energy = float(
        0.5
        * np.mean(
            np.fft.ifftn(legal_u0_hat[0]).real ** 2
            + np.fft.ifftn(legal_u0_hat[1]).real ** 2
            + np.fft.ifftn(legal_u0_hat[2]).real ** 2
        )
        * (2 * np.pi) ** 3
    )
    illegal_u0_hat = rescale_energy(illegal_shell, target_energy=legal_energy)

    legal_cert = _cert_from_u0(
        legal_u0_hat,
        N=args.N,
        window_c=args.window_c,
        energy_threshold=args.energy_threshold,
        kappa_threshold=args.kappa_threshold,
        burst_threshold=args.burst_threshold,
        config={"case": "LEGAL", **cfg},
    )
    illegal_cert = _cert_from_u0(
        illegal_u0_hat,
        N=args.N,
        window_c=args.window_c,
        energy_threshold=args.energy_threshold,
        kappa_threshold=args.kappa_threshold,
        burst_threshold=args.burst_threshold,
        config={"case": "ILLEGAL", **cfg},
    )

    legal_max_kappa = float(legal_cert["summary"]["max_kappa_max_over_time"] or 0.0)
    legal_max_burst = float(legal_cert["summary"]["max_burst_ratio"] or 0.0)
    illegal_max_kappa = float(illegal_cert["summary"]["max_kappa_max_over_time"] or 0.0)
    illegal_max_burst = float(illegal_cert["summary"]["max_burst_ratio"] or 0.0)

    print(
        f"CASE=LEGAL verdict={legal_cert['verdict']} "
        f"max_kappa={legal_max_kappa:.6e} max_burst={legal_max_burst:.6e}"
    )
    print(
        f"CASE=ILLEGAL verdict={illegal_cert['verdict']} "
        f"max_kappa={illegal_max_kappa:.6e} max_burst={illegal_max_burst:.6e}"
    )

    j_legal = np.asarray(legal_cert["metrics"]["j_values"], dtype=int)
    j_illegal = np.asarray(illegal_cert["metrics"]["j_values"], dtype=int)
    kappa_legal = np.asarray(legal_cert["metrics"]["kappa_max_over_time"], dtype=float)
    kappa_illegal = np.asarray(illegal_cert["metrics"]["kappa_max_over_time"], dtype=float)

    fig = plt.figure(figsize=(6.4, 4.2))
    plt.plot(j_legal, kappa_legal, "o-", label="LEGAL")
    plt.plot(j_illegal, kappa_illegal, "s-", label="ILLEGAL")
    plt.axhline(args.kappa_threshold, color="r", linestyle="--", linewidth=1.0, label="kappa_threshold")
    plt.xlabel("dyadic shell j")
    plt.ylabel("kappa_max_over_time(j) at t=0")
    plt.title("SPT legal vs intentionally illegal initial condition")
    plt.grid(alpha=0.3)
    plt.legend()
    fig.tight_layout()

    png_path = out_dir / "ns_spt_illegal_vs_legal.png"
    plt.savefig(png_path, dpi=150)
    print(f"saved {png_path}")

    npz_path = out_dir / "ns_spt_illegal_vs_legal_data.npz"
    meta_json = np.frombuffer(json.dumps(cfg, sort_keys=True).encode("utf-8"), dtype=np.uint8)
    savez_deterministic(
        npz_path,
        meta_json=meta_json,
        j_legal=j_legal,
        j_illegal=j_illegal,
        kappa_legal=kappa_legal,
        kappa_illegal=kappa_illegal,
        burst_legal=np.asarray(legal_cert["metrics"]["burst_ratio"], dtype=float),
        burst_illegal=np.asarray(illegal_cert["metrics"]["burst_ratio"], dtype=float),
        legal_verdict=np.array(legal_cert["verdict"]),
        illegal_verdict=np.array(illegal_cert["verdict"]),
        legal_max_kappa=np.array(legal_max_kappa),
        illegal_max_kappa=np.array(illegal_max_kappa),
        legal_max_burst=np.array(legal_max_burst),
        illegal_max_burst=np.array(illegal_max_burst),
        legal_threshold=np.array(args.kappa_threshold),
    )
    print(f"saved {npz_path}")

    expected_ok = (legal_cert["verdict"] == "PASS") and (illegal_cert["verdict"] == "FAIL")
    return 0 if expected_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
