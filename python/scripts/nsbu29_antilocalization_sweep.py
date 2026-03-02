import argparse
import json
import sys
from dataclasses import dataclass
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
    u_shell = np.empty_like(u_hat.real)
    for i in range(3):
        u_shell[i] = np.fft.ifftn(u_shell_hat[i]).real
    return u_shell


def m95_from_snapshot_matrix(X: np.ndarray, energy_threshold: float) -> int:
    x_norm = float(np.linalg.norm(X))
    if x_norm <= 1e-14:
        return 0
    svals = np.linalg.svd(X, full_matrices=False, compute_uv=False)
    s2 = svals * svals
    total = float(np.sum(s2))
    if total <= 0.0:
        return 0
    cum = np.cumsum(s2) / total
    return int(np.searchsorted(cum, energy_threshold) + 1)


@dataclass
class CaseSummary:
    j_values: np.ndarray
    kappa_max: np.ndarray
    burst_ratio: np.ndarray
    m95: np.ndarray


def run_case(
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
) -> CaseSummary:
    L = 2 * np.pi
    u0_hat = make_deterministic_u0_hat_3d(N, L=L, seed=seed, k_init=k_init, dealias=True)
    f_hat = make_lowk_divfree_forcing_hat_3d(
        N,
        L=L,
        seed=seed,
        k_force=k_force,
        amp=forcing_amp,
        dealias=True,
    )

    def forcing_hat_fn(_t: float, _kx: np.ndarray, _ky: np.ndarray, _kz: np.ndarray) -> np.ndarray:
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
    if res.u_hat_hist is None or res.t_hat is None:
        raise RuntimeError("u_hat_hist/t_hat were not recorded")

    _, _, _, k2 = kgrid_3d(N, L=L)
    k_mag = np.sqrt(k2)
    k_max = float(np.max(k_mag))
    j_max = int(np.floor(np.log2(max(k_max, 1.0))))
    j_values = np.arange(1, j_max + 1, dtype=int)

    shell_masks = [dyadic_shell_mask(k_mag, int(j), k0=1.0) for j in j_values]
    windows = [window_for_shell(N, int(j), c=window_c) for j in j_values]

    n_snap = res.u_hat_hist.shape[0]
    n_feat = 3 * N * N * N
    conc = np.zeros((n_snap, len(j_values)), dtype=float)
    m95 = np.zeros(len(j_values), dtype=int)

    for jj, (mask, window) in enumerate(zip(shell_masks, windows)):
        X = np.empty((n_snap, n_feat), dtype=np.float32)
        for i in range(n_snap):
            u_shell = shell_real_space(res.u_hat_hist[i], mask)
            conc[i, jj] = max_local_energy_fraction(u_shell, window)
            X[i] = u_shell.reshape(-1)
        m95[jj] = m95_from_snapshot_matrix(X, energy_threshold=energy_threshold)

    mean_conc = np.mean(conc, axis=0)
    max_conc = np.max(conc, axis=0)
    _, _, kappa_mean = localization_factors(
        mean_conc,
        N=N,
        windows=np.array(windows),
        j_values=j_values,
    )
    _, _, kappa_max = localization_factors(
        max_conc,
        N=N,
        windows=np.array(windows),
        j_values=j_values,
    )
    burst_ratio = max_conc / np.maximum(mean_conc, 1e-30)

    _ = kappa_mean  # retained for local debug symmetry
    return CaseSummary(
        j_values=j_values,
        kappa_max=kappa_max,
        burst_ratio=burst_ratio,
        m95=m95,
    )


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="NS-BU-29 anti-localization stress sweep")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--resolutions", type=int, nargs="+", default=[32, 48])
    p.add_argument("--nu-list", type=float, nargs="+", default=[0.004, 0.002, 0.001])
    p.add_argument("--seeds", type=int, nargs="+", default=None)
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
    p.add_argument("--out-dir", default=str(ROOT / "python" / "artifacts" / "paper"))
    return p.parse_args()


def main() -> int:
    args = parse_args()
    seeds = args.seeds if args.seeds is not None else [args.seed, args.seed + 1, args.seed + 2]
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cfg = {
        "script": Path(__file__).name,
        "seed": args.seed,
        "resolutions": args.resolutions,
        "nu_list": args.nu_list,
        "seeds": seeds,
        "dt": args.dt,
        "t_max": args.t_max,
        "nu": None,
        "mu": args.mu,
        "alpha": args.alpha,
        "k_init": args.k_init,
        "k_force": args.k_force,
        "amp": args.forcing_amp,
        "record_u_hat_every": args.record_u_hat_every,
        "window_c": args.window_c,
        "energy_threshold": args.energy_threshold,
        "kappa_threshold": args.kappa_threshold,
        "burst_threshold": args.burst_threshold,
        "out_dir": str(out_dir),
    }
    print(f"FIGURE_CONFIG: {json.dumps(cfg, sort_keys=True)}")

    combos = [(N, nu, seed) for N in args.resolutions for nu in args.nu_list for seed in seeds]
    total_runs = len(combos)
    print("NS-BU-29 anti-localization stress sweep")
    print(
        f"runs={total_runs} resolutions={args.resolutions} nu_list={args.nu_list} "
        f"seeds={seeds} dt={args.dt} t_max={args.t_max}"
    )

    agg: dict[int, dict[str, list[float]]] = {}
    failures: list[tuple[int, float, int, str]] = []

    for idx, (N, nu, seed) in enumerate(combos, start=1):
        print(f"\n[{idx}/{total_runs}] N={N} nu={nu:.4g} seed={seed}")
        try:
            case = run_case(
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
        except Exception as exc:  # pragma: no cover - sweep resilience path
            failures.append((N, nu, seed, str(exc)))
            print(f"  FAIL: {exc}")
            continue

        for j, kappa, burst, m95 in zip(
            case.j_values,
            case.kappa_max,
            case.burst_ratio,
            case.m95,
        ):
            j_int = int(j)
            if j_int not in agg:
                agg[j_int] = {"kappa": [], "burst": [], "m95": []}
            agg[j_int]["kappa"].append(float(kappa))
            agg[j_int]["burst"].append(float(burst))
            agg[j_int]["m95"].append(float(m95))

        print(
            f"  done: max_kappa={float(np.max(case.kappa_max)):.6e} "
            f"max_burst={float(np.max(case.burst_ratio)):.6e} "
            f"max_m95={int(np.max(case.m95))}"
        )

    if not agg:
        print("No successful runs; aborting.")
        return 1

    j_sorted = np.array(sorted(agg.keys()), dtype=int)
    worst_kappa = np.array([max(agg[int(j)]["kappa"]) for j in j_sorted], dtype=float)
    worst_burst = np.array([max(agg[int(j)]["burst"]) for j in j_sorted], dtype=float)
    worst_m95 = np.array([max(agg[int(j)]["m95"]) for j in j_sorted], dtype=float)
    n_cov = np.array([len(agg[int(j)]["kappa"]) for j in j_sorted], dtype=int)

    print("\nWorst-case summary by shell:")
    print("j | runs | worst_kappa_max_over_time | worst_burst_ratio | worst_m95")
    for j, n, k, b, m in zip(j_sorted, n_cov, worst_kappa, worst_burst, worst_m95):
        print(f"{j:2d} | {n:4d} | {k:25.6e} | {b:17.6e} | {int(m):9d}")

    max_kappa = float(np.max(worst_kappa))
    max_burst = float(np.max(worst_burst))
    verdict = (max_kappa < args.kappa_threshold) and (max_burst < args.burst_threshold)
    print(
        "ANTI_LOCALIZATION_PLAUSIBLE = "
        f"{verdict} (max_kappa={max_kappa:.6e} < {args.kappa_threshold}, "
        f"max_burst={max_burst:.6e} < {args.burst_threshold})"
    )

    if failures:
        print("\nFailed runs:")
        for N, nu, seed, msg in failures:
            print(f"  N={N} nu={nu} seed={seed}: {msg}")

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    axes[0].plot(j_sorted, worst_kappa, "o-")
    axes[0].axhline(args.kappa_threshold, color="r", linestyle="--", linewidth=1.0)
    axes[0].set_xlabel("dyadic shell j")
    axes[0].set_ylabel("worst kappa_max_over_time")
    axes[0].set_title("Worst Kappa by Shell")
    axes[0].grid(alpha=0.3)

    axes[1].plot(j_sorted, worst_burst, "s-")
    axes[1].axhline(args.burst_threshold, color="r", linestyle="--", linewidth=1.0)
    axes[1].set_xlabel("dyadic shell j")
    axes[1].set_ylabel("worst burst_ratio")
    axes[1].set_title("Worst Burst Ratio by Shell")
    axes[1].grid(alpha=0.3)

    fig.suptitle("NS-BU-29 anti-localization stress sweep")
    fig.tight_layout()
    out_path = out_dir / "nsbu29_antilocalization_sweep.png"
    plt.savefig(out_path, dpi=150)
    print(f"saved {out_path}")

    data_path = out_dir / "nsbu29_antilocalization_sweep_data.npz"
    meta_json = np.frombuffer(json.dumps(cfg, sort_keys=True).encode("utf-8"), dtype=np.uint8)
    failure_arr = np.array(
        [f"N={N} nu={nu} seed={seed}: {msg}" for (N, nu, seed, msg) in failures],
        dtype=str,
    )
    savez_deterministic(
        data_path,
        meta_json=meta_json,
        j_sorted=j_sorted,
        worst_kappa=worst_kappa,
        worst_burst=worst_burst,
        worst_m95=worst_m95,
        n_cov=n_cov,
        max_kappa=np.array(max_kappa),
        max_burst=np.array(max_burst),
        verdict=np.array(int(verdict)),
        failures=failure_arr,
    )
    print(f"saved data: {data_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
