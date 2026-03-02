import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "python" / "src"
sys.path.insert(0, str(SRC))

from nswave.dyadic import dyadic_shell_mask  # noqa: E402
from nswave.localization import (  # noqa: E402
    localization_factors,
    max_local_energy_fraction,
    window_for_shell,
)
from nswave.ns3d_spectral import make_deterministic_u0_hat_3d, simulate_ns3d  # noqa: E402
from nswave.operators import dealias_mask_3d, kgrid_3d, project_div_free_3d  # noqa: E402
from nswave.deterministic_npz import savez_deterministic  # noqa: E402

ART = ROOT / "python" / "artifacts"


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


def empirical_channelization_summary(
    u_hat_hist: np.ndarray,
    shell_masks: list[np.ndarray],
    windows: list[int],
    energy_threshold: float,
    top_modes: int,
) -> list[tuple[int, int, float]]:
    n_snap = u_hat_hist.shape[0]
    N = u_hat_hist.shape[-1]
    n_feat = 3 * N * N * N
    rows: list[tuple[int, int, float]] = []

    print("\nEmpirical channelization check (diagnostic-only):")
    for j, (mask, window) in enumerate(zip(shell_masks, windows)):
        X = np.empty((n_snap, n_feat), dtype=float)
        for i in range(n_snap):
            u_shell = shell_real_space(u_hat_hist[i], mask)
            X[i] = u_shell.reshape(-1)

        X_norm = float(np.linalg.norm(X))
        if X_norm <= 1e-14:
            print(f"  j={j:2d} m95=0 (inactive shell) worst_top_loc=nan")
            rows.append((j, 0, float("nan")))
            continue

        _, svals, vh = np.linalg.svd(X, full_matrices=False)
        s2 = svals * svals
        total = float(np.sum(s2))
        if total <= 0.0:
            print(f"  j={j:2d} m95=0 (degenerate shell) worst_top_loc=nan")
            rows.append((j, 0, float("nan")))
            continue

        cum = np.cumsum(s2) / total
        m95 = int(np.searchsorted(cum, energy_threshold) + 1)

        top_k = min(top_modes, vh.shape[0])
        loc_vals = []
        for r in range(top_k):
            mode = vh[r].reshape(3, N, N, N)
            loc_vals.append(max_local_energy_fraction(mode, window))
        worst_top_loc = float(np.max(loc_vals)) if loc_vals else float("nan")
        print(
            f"  j={j:2d} m95={m95:3d}/{n_snap:3d} "
            f"worst_top_loc={worst_top_loc:.6e}"
        )
        rows.append((j, m95, worst_top_loc))

    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="NS-BU-28 anti-localization diagnostic")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--N", type=int, default=32)
    parser.add_argument("--dt", type=float, default=5e-4)
    parser.add_argument("--t-max", type=float, default=0.6)
    parser.add_argument("--nu", type=float, default=1e-3)
    parser.add_argument("--mu", type=float, default=0.0)
    parser.add_argument("--alpha", type=float, default=2.0)
    parser.add_argument("--k-init", type=int, default=10)
    parser.add_argument("--k-force", type=int, default=2)
    parser.add_argument("--forcing-amp", type=float, default=1e-2)
    parser.add_argument("--seed-u0", type=int, default=None)
    parser.add_argument("--seed-force", type=int, default=None)
    parser.add_argument("--record-u-hat-every", type=int, default=20)
    parser.add_argument("--window-c", type=float, default=1.0)
    parser.add_argument("--energy-threshold", type=float, default=0.95)
    parser.add_argument("--top-modes", type=int, default=3)
    parser.add_argument("--out-dir", default=str(ROOT / "python" / "artifacts" / "paper"))
    parser.add_argument("--skip-channelization", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    L = 2 * np.pi
    seed_u0 = args.seed if args.seed_u0 is None else args.seed_u0
    seed_force = args.seed if args.seed_force is None else args.seed_force
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cfg = {
        "script": Path(__file__).name,
        "seed": args.seed,
        "seed_u0": seed_u0,
        "seed_force": seed_force,
        "N": args.N,
        "dt": args.dt,
        "t_max": args.t_max,
        "nu": args.nu,
        "mu": args.mu,
        "alpha": args.alpha,
        "k_init": args.k_init,
        "k_force": args.k_force,
        "amp": args.forcing_amp,
        "record_u_hat_every": args.record_u_hat_every,
        "window_c": args.window_c,
        "energy_threshold": args.energy_threshold,
        "top_modes": args.top_modes,
        "skip_channelization": bool(args.skip_channelization),
        "out_dir": str(out_dir),
    }
    print(f"FIGURE_CONFIG: {json.dumps(cfg, sort_keys=True)}")

    u0_hat = make_deterministic_u0_hat_3d(
        args.N,
        L=L,
        seed=seed_u0,
        k_init=args.k_init,
        dealias=True,
    )
    f_hat = make_lowk_divfree_forcing_hat_3d(
        args.N,
        L=L,
        seed=seed_force,
        k_force=args.k_force,
        amp=args.forcing_amp,
        dealias=True,
    )

    def forcing_hat_fn(_t: float, _kx: np.ndarray, _ky: np.ndarray, _kz: np.ndarray) -> np.ndarray:
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
    if res.u_hat_hist is None or res.t_hat is None:
        raise RuntimeError("u_hat_hist/t_hat were not recorded")

    _, _, _, k2 = kgrid_3d(args.N, L=L)
    k_mag = np.sqrt(k2)
    k_max = float(np.max(k_mag))
    j_max = int(np.floor(np.log2(max(k_max, 1.0))))
    j_values = np.arange(0, j_max + 1, dtype=int)

    shell_masks = [dyadic_shell_mask(k_mag, int(j), k0=1.0) for j in j_values]
    windows = [window_for_shell(args.N, int(j), c=args.window_c) for j in j_values]

    n_snap = res.u_hat_hist.shape[0]
    conc = np.zeros((n_snap, len(j_values)), dtype=float)
    for i in range(n_snap):
        u_hat = res.u_hat_hist[i]
        for jj, (mask, window) in enumerate(zip(shell_masks, windows)):
            u_shell = shell_real_space(u_hat, mask)
            conc[i, jj] = max_local_energy_fraction(u_shell, window)

    mean_conc = np.mean(conc, axis=0)
    max_conc = np.max(conc, axis=0)
    f_j, L_mean, kappa_mean = localization_factors(
        mean_conc,
        N=args.N,
        windows=np.array(windows),
        j_values=j_values,
    )
    _, L_max, kappa_max = localization_factors(
        max_conc,
        N=args.N,
        windows=np.array(windows),
        j_values=j_values,
    )

    print("NS-BU-28 anti-localization diagnostic")
    print(
        f"N={args.N} dt={args.dt} t_max={args.t_max} nu={args.nu} "
        f"mu={args.mu} alpha={args.alpha} k_force={args.k_force} amp={args.forcing_amp}"
    )
    print(f"snapshots={n_snap} record_u_hat_every={args.record_u_hat_every}")
    print("j | window | mean_conc | max_conc | L_mean | L_max | kappa_mean | kappa_max")
    for j, w, mc, xc, lm, lx, km, kx in zip(
        j_values,
        windows,
        mean_conc,
        max_conc,
        L_mean,
        L_max,
        kappa_mean,
        kappa_max,
    ):
        print(
            f"{j:2d} | {w:6d} | {mc:9.6e} | {xc:9.6e} | "
            f"{lm:9.6e} | {lx:9.6e} | {km:9.6e} | {kx:9.6e}"
        )

    if len(j_values) >= 2:
        diffs = np.diff(mean_conc)
        monotone_dec = bool(np.all(diffs <= 1e-4))
        slope = float(np.polyfit(j_values, mean_conc, 1)[0])
    else:
        monotone_dec = True
        slope = 0.0
    print(
        "summary: mean_conc decreases with j (heuristic) = "
        f"{monotone_dec} (linear slope={slope:.6e})"
    )

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))

    axes[0].plot(j_values, mean_conc, "o-", label="measured mean concentration")
    ref = (j_values + 1.0) * (2.0 ** (-3.0 * j_values))
    if ref[0] > 0:
        ref = ref * (mean_conc[0] / ref[0])
    axes[0].plot(j_values, ref, "--", label="heuristic ~(j+1)2^{-3j} (rescaled)")
    axes[0].set_yscale("log")
    axes[0].set_xlabel("dyadic shell j")
    axes[0].set_ylabel("mean max-local-energy fraction")
    axes[0].set_title("Concentration")
    axes[0].grid(alpha=0.3)
    axes[0].legend()

    axes[1].plot(j_values, kappa_mean, "o-", label="kappa_mean")
    axes[1].plot(j_values, kappa_max, "s--", label="kappa_max")
    axes[1].set_xlabel("dyadic shell j")
    axes[1].set_ylabel("kappa")
    axes[1].set_title("Channel-Compatible Kappa")
    axes[1].grid(alpha=0.3)
    axes[1].legend()

    fig.suptitle("NS-BU-28 anti-localization diagnostic")
    fig.tight_layout()

    out_path = out_dir / "nsbu28_antilocalization.png"
    plt.savefig(out_path, dpi=150)
    print(f"saved {out_path}")

    channel_rows: list[tuple[int, int, float]] = []
    if not args.skip_channelization:
        channel_rows = empirical_channelization_summary(
            res.u_hat_hist,
            shell_masks,
            windows,
            energy_threshold=args.energy_threshold,
            top_modes=args.top_modes,
        )

    data_path = out_dir / "nsbu28_antilocalization_data.npz"
    meta_json = np.frombuffer(json.dumps(cfg, sort_keys=True).encode("utf-8"), dtype=np.uint8)
    channel_arr = (
        np.array(channel_rows, dtype=float) if channel_rows else np.empty((0, 3), dtype=float)
    )
    savez_deterministic(
        data_path,
        meta_json=meta_json,
        j_values=j_values,
        windows=np.array(windows, dtype=int),
        conc=conc,
        mean_conc=mean_conc,
        max_conc=max_conc,
        f_j=f_j,
        L_mean=L_mean,
        L_max=L_max,
        kappa_mean=kappa_mean,
        kappa_max=kappa_max,
        channel_rows=channel_arr,
    )
    print(f"saved data: {data_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
