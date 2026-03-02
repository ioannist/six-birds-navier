import argparse
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "python" / "src"
sys.path.insert(0, str(SRC))

from nswave.capacity_estimator import capacity_L2_ratio
from nswave.coarse_grain import restrict_hat_3d
from nswave.dyadic import dyadic_shell_mask
from nswave.deterministic_npz import savez_deterministic
from nswave.ns3d_spectral import make_deterministic_u0_hat_3d, nonlinear_term_hat_3d, simulate_ns3d
from nswave.operators import dealias_mask_3d, kgrid_3d, project_div_free_3d


ART = ROOT / "python" / "artifacts"


def make_lowk_divfree_forcing_hat_3d(
    N: int,
    *,
    L: float,
    seed: int,
    k_force: int,
    amp: float,
    dealias: bool,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    f_phys = rng.standard_normal((3, N, N, N))
    f_hat = np.zeros_like(f_phys, dtype=complex)
    for i in range(3):
        f_hat[i] = np.fft.fftn(f_phys[i])
    f_hat[:, 0, 0, 0] = 0.0

    kk = np.fft.fftfreq(N) * N
    kx_int, ky_int, kz_int = np.meshgrid(kk, kk, kk, indexing="ij")
    k_mag_int = np.sqrt(kx_int**2 + ky_int**2 + kz_int**2)
    f_hat *= k_mag_int <= k_force

    kx, ky, kz, _ = kgrid_3d(N, L=L)
    f_hat = project_div_free_3d(f_hat, kx, ky, kz)
    if dealias:
        f_hat *= dealias_mask_3d(N)

    f_phys = [np.fft.ifftn(f_hat[i]).real for i in range(3)]
    rms = np.sqrt(np.mean(f_phys[0] ** 2 + f_phys[1] ** 2 + f_phys[2] ** 2))
    scale = amp / max(rms, 1e-12)
    return f_hat * scale


def main() -> int:
    parser = argparse.ArgumentParser(description="NS-BU-21 route-capacity numbers")
    parser.add_argument("--quick", action="store_true", help="Run a short quick configuration")
    parser.add_argument("--ladder", default="default", choices=["default", "auto_kcut"], help="Ladder mode")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out-dir", default=str(ROOT / "python" / "artifacts" / "paper"))
    args = parser.parse_args()

    L = 2 * np.pi
    dt = 5e-4
    nu = 1e-3
    mu = 0.0
    alpha = 2.0
    dealias = True
    seed = args.seed
    k_init = 3
    record_u_hat_every = 20
    burn_in_frac = 0.4
    shell_j = 1
    k_force = 2
    amp = 1.0e-2
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.quick:
        N_levels = [8, 16, 32]
        t_max = 0.05
        record_u_hat_every = 10
    else:
        N_levels = [8, 12, 16, 24, 32]
        t_max = 0.2

    if args.ladder == "auto_kcut":
        N0 = 8
        N_hi = max(48, N_levels[-1])
        auto_levels = [N0]
        kx0, ky0, kz0, k2_0 = kgrid_3d(N0, L=L)
        k_mag0 = np.sqrt(k2_0)
        mask0 = dealias_mask_3d(N0).astype(bool) if dealias else np.ones_like(k_mag0, dtype=bool)
        prev_cut = float(np.max(k_mag0[mask0]))
        for N in range(N0 + 2, N_hi + 1, 2):
            kx, ky, kz, k2 = kgrid_3d(N, L=L)
            k_mag = np.sqrt(k2)
            mask = dealias_mask_3d(N).astype(bool) if dealias else np.ones_like(k_mag, dtype=bool)
            k_cut = float(np.max(k_mag[mask]))
            if k_cut > prev_cut + 1e-12:
                auto_levels.append(N)
                prev_cut = k_cut
        if auto_levels[-1] != N_hi:
            auto_levels.append(N_hi)
        N_levels = auto_levels

    N_hi = N_levels[-1]
    N0 = N_levels[0]
    cfg = {
        "script": Path(__file__).name,
        "seed": seed,
        "N": N_hi,
        "dt": dt,
        "t_max": t_max,
        "nu": nu,
        "mu": mu,
        "alpha": alpha,
        "N_levels": N_levels,
        "shell_j": shell_j,
        "burn_in_frac": burn_in_frac,
        "k_force": k_force,
        "amp": amp,
        "k_init": k_init,
        "out_dir": str(out_dir),
    }
    print(f"FIGURE_CONFIG: {json.dumps(cfg, sort_keys=True)}")

    u0_hat = make_deterministic_u0_hat_3d(
        N_hi, L=L, seed=seed, k_init=k_init, dealias=dealias
    )
    forcing_hat_hi = make_lowk_divfree_forcing_hat_3d(
        N_hi, L=L, seed=seed, k_force=k_force, amp=amp, dealias=dealias
    )

    def forcing_hat_fn(t: float, kx: np.ndarray, ky: np.ndarray, kz: np.ndarray) -> np.ndarray:
        return forcing_hat_hi

    res = simulate_ns3d(
        N=N_hi,
        dt=dt,
        t_max=t_max,
        nu=nu,
        mu=mu,
        alpha=alpha,
        u0_hat=u0_hat,
        forcing_hat_fn=forcing_hat_fn,
        dealias=dealias,
        record_every=10,
        record_u_hat_every=record_u_hat_every,
    )

    if res.u_hat_hist is None or res.t_hat is None:
        raise RuntimeError("u_hat_hist not recorded")

    burn_idx = int(burn_in_frac * len(res.t_hat))
    u_hist_hi = res.u_hat_hist[burn_idx:]

    dt_hat = dt * record_u_hat_every

    kx0, ky0, kz0, k2_0 = kgrid_3d(N0, L=L)
    mask0 = dealias_mask_3d(N0) if dealias else None
    k_mag0 = np.sqrt(k2_0)
    shell_mask = dyadic_shell_mask(k_mag0, shell_j, k0=1.0)

    u0_hist = np.empty((len(u_hist_hi), 3, N0, N0, N0), dtype=complex)
    for idx, u_hi in enumerate(u_hist_hi):
        u0 = restrict_hat_3d(u_hi, N0)
        if mask0 is not None:
            u0 *= mask0
        u0_hist[idx] = u0

    u_rms_samples = []
    for u_hat in u0_hist:
        u = [np.fft.ifftn(u_hat[i]).real for i in range(3)]
        mean_u2 = np.mean(u[0] ** 2 + u[1] ** 2 + u[2] ** 2)
        u_rms_samples.append(mean_u2)
    u_rms = float(np.sqrt(np.mean(u_rms_samples)))

    nl0_hist = []
    for N_level in N_levels:
        kx, ky, kz, k2 = kgrid_3d(N_level, L=L)
        mask = dealias_mask_3d(N_level) if dealias else None
        nl0 = np.empty((len(u_hist_hi), 3, N0, N0, N0), dtype=complex)
        for idx, u_hi in enumerate(u_hist_hi):
            u_lvl = restrict_hat_3d(u_hi, N_level)
            if mask is not None:
                u_lvl *= mask
            nl_lvl = nonlinear_term_hat_3d(u_lvl, kx, ky, kz, k2, mask=mask)
            nl0[idx] = restrict_hat_3d(nl_lvl, N0)
        nl0_hist.append(nl0)

    u_hist = u0_hist * shell_mask
    u_flat = u_hist.reshape(u_hist.shape[0], -1)
    u_port_energy = float(np.sum(np.abs(u_flat) ** 2) * dt_hat)
    n_levels = len(N_levels)
    n_increments = n_levels - 1
    lambda_inc = np.zeros(n_increments)
    y_inc_list = []

    for m in range(n_increments):
        c_inc = -(nl0_hist[m + 1] - nl0_hist[m])
        y_inc = c_inc * shell_mask
        y_inc_list.append(y_inc)
        lambda_inc[m] = capacity_L2_ratio(u_hist, y_inc, dt_hat)

    c_sum = np.zeros(n_levels)
    for n in range(1, n_levels):
        c_sum[n] = np.sum(lambda_inc[:n])

    lambda_direct = np.zeros(n_levels)
    e_route = np.zeros(n_levels)
    for n in range(1, n_levels):
        c_direct = -(nl0_hist[n] - nl0_hist[0])
        y_direct = c_direct * shell_mask
        lambda_direct[n] = capacity_L2_ratio(u_hist, y_direct, dt_hat)
        y_sum = np.zeros_like(y_direct)
        for m in range(n):
            y_sum += y_inc_list[m]
        y_route = y_direct - y_sum
        e_route[n] = capacity_L2_ratio(u_hist, y_route, dt_hat)

    eps = 1e-12
    e_vals = np.abs(lambda_direct - c_sum)
    slack_poswork = np.maximum(0.0, c_sum - lambda_direct)
    partial_sum_e = np.cumsum(e_vals)
    partial_sum_e_route = np.cumsum(e_route)
    r_vals = np.full(n_levels, np.nan)
    for n in range(1, n_levels - 1):
        r_vals[n] = c_sum[n + 1] / max(c_sum[n], eps)

    k_cut = np.zeros(n_levels)
    for idx, N_level in enumerate(N_levels):
        kx, ky, kz, k2 = kgrid_3d(N_level, L=L)
        k_mag = np.sqrt(k2)
        if dealias:
            mask = dealias_mask_3d(N_level).astype(bool)
        else:
            mask = np.ones_like(k_mag, dtype=bool)
        k_cut[idx] = float(np.max(k_mag[mask]))

    c_tilde = np.zeros(n_levels)
    for n in range(1, n_levels):
        denom = max(u_rms * k_cut[n], eps)
        c_tilde[n] = c_sum[n] / denom

    r_tilde = np.full(n_levels, np.nan)
    for n in range(2, n_levels):
        r_tilde[n] = c_tilde[n] / max(c_tilde[n - 1], eps)

    p_candidates = [0, 1, 2, 3]
    p_fit = None
    for p in p_candidates:
        ok = True
        for n in range(1, n_levels):
            if not np.isfinite(r_tilde[n]):
                continue
            r_model = ((n + 1) / n) ** p
            if r_tilde[n] > r_model + 1e-12:
                ok = False
                break
        if ok:
            p_fit = p
            break

    print("NS-BU-21 route-capacity numbers")
    print(f"N_hi={N_hi} N_levels={N_levels} shell_j={shell_j} burn_in_frac={burn_in_frac}")
    print(f"n_steps={n_levels - 1}")
    print(f"k_cut={k_cut.tolist()}")
    print(f"U_rms={u_rms:.6e}")
    print("n | N_level | k_cut | Lambda_direct | C_sum | C_tilde | r_tilde | e_route | slack_poswork")
    for n in range(1, n_levels):
        r_val = r_tilde[n]
        r_txt = f"{r_val:.6e}" if np.isfinite(r_val) else "-"
        print(
            f"{n} | {N_levels[n]} | {k_cut[n]:.6e} | {lambda_direct[n]:.6e} | {c_sum[n]:.6e} | "
            f"{c_tilde[n]:.6e} | {r_txt} | {e_route[n]:.6e} | {slack_poswork[n]:.6e}"
        )

    if np.max(e_route[1:]) <= 1e-12:
        print(
            "NOTE: route mismatch is ~0 (telescoping closure + associative restriction); "
            "mismatch-summability is therefore trivially satisfied for this ladder."
        )

    print(
        "NOTE: C_direct telescopes exactly at the field level; e(n) measures capacity-currency "
        "non-additivity/slack, not algebraic mismatch."
    )
    print(f"p_fit={p_fit}")
    print(f"max_e_route={float(np.max(e_route[1:])):.6e}")
    print(
        f"slack_poswork_min={float(np.min(slack_poswork[1:])):.6e} "
        f"slack_poswork_max={float(np.max(slack_poswork[1:])):.6e}"
    )

    fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=False)
    n_vals = np.arange(1, n_levels)
    axes[0, 0].plot(n_vals, np.maximum(c_sum[1:], eps), marker="o", label="C_sum")
    axes[0, 0].plot(n_vals, np.maximum(lambda_direct[1:], eps), marker="s", label="Lambda_direct")
    axes[0, 0].set_yscale("log")
    axes[0, 0].set_ylabel("capacity (log)")
    axes[0, 0].legend(loc="best")

    axes[1, 0].plot(n_vals, e_vals[1:], marker="o", label="e(n)")
    axes[1, 0].plot(n_vals, partial_sum_e[1:], marker="s", label="sum e")
    axes[1, 0].set_xlabel("depth n")
    axes[1, 0].set_ylabel("mismatch")
    axes[1, 0].legend(loc="best")

    axes[0, 1].plot(n_vals, np.maximum(c_tilde[1:], eps), marker="o", label="C_tilde")
    axes[0, 1].set_yscale("log")
    axes[0, 1].set_ylabel("C_tilde (log)")
    axes[0, 1].legend(loc="best")

    axes[1, 1].plot(n_vals, r_tilde[1:], marker="o", label="r_tilde")
    for p in p_candidates:
        model = [((n + 1) / n) ** p for n in n_vals]
        axes[1, 1].plot(n_vals, model, linestyle="--", label=f"p={p}")
    axes[1, 1].set_xlabel("depth n")
    axes[1, 1].set_ylabel("r_tilde")
    axes[1, 1].legend(loc="best")

    plt.tight_layout()
    out_path = out_dir / "nsbu21_route_capacity_numbers.png"
    plt.savefig(out_path, dpi=150)
    print(f"saved {out_path}")

    e_tilde = e_vals / np.maximum(u_rms * k_cut, eps)
    e_route_tilde = e_route / np.maximum(u_rms * k_cut, eps)
    data_path = out_dir / "nsbu21_route_capacity_numbers_data.npz"
    meta_json = np.frombuffer(json.dumps(cfg, sort_keys=True).encode("utf-8"), dtype=np.uint8)
    savez_deterministic(
        data_path,
        meta_json=meta_json,
        n=np.arange(n_levels, dtype=int),
        N_levels=np.array(N_levels, dtype=int),
        k_cut=k_cut,
        Lambda_direct=lambda_direct,
        C_sum=c_sum,
        e=e_vals,
        e_tilde=e_tilde,
        e_route=e_route,
        e_route_tilde=e_route_tilde,
        slack_poswork=slack_poswork,
        partial_sum_e=partial_sum_e,
        partial_sum_e_route=partial_sum_e_route,
        U_rms=np.array(u_rms),
        C_tilde=c_tilde,
        r_tilde=r_tilde,
        u_port_energy=np.full(n_levels, u_port_energy),
        N_hi=np.array(N_hi),
        dt=np.array(dt),
        t_max=np.array(t_max),
        nu=np.array(nu),
        mu=np.array(mu),
        alpha=np.array(alpha),
        shell_j=np.array(shell_j),
        burn_in_frac=np.array(burn_in_frac),
        k_force=np.array(k_force),
        forcing_amp=np.array(amp),
        seed=np.array(seed),
        k_init=np.array(k_init),
    )
    print(f"saved data: {data_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
