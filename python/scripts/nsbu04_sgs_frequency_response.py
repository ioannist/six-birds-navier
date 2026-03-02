import argparse
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "python" / "src"
sys.path.insert(0, str(SRC))

from nswave.coarse_grain import restrict_hat_3d  # noqa: E402
from nswave.deterministic_npz import savez_deterministic  # noqa: E402
from nswave.dyadic import dyadic_shell_index  # noqa: E402
from nswave.ns3d_spectral import make_deterministic_u0_hat_3d, simulate_ns3d, nonlinear_term_hat_3d  # noqa: E402
from nswave.operators import dealias_mask_3d, kgrid_3d, project_div_free_3d  # noqa: E402


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

    kx_int = np.fft.fftfreq(N) * N
    ky_int = kx_int.copy()
    kz_int = kx_int.copy()
    kx, ky, kz = np.meshgrid(kx_int, ky_int, kz_int, indexing="ij")
    k_mag_int = np.sqrt(kx**2 + ky**2 + kz**2)
    mask = k_mag_int <= k_force
    f_hat *= mask

    kx_phys, ky_phys, kz_phys, k2 = kgrid_3d(N, L=L)
    f_hat = project_div_free_3d(f_hat, kx_phys, ky_phys, kz_phys)
    if dealias:
        f_hat *= dealias_mask_3d(N)

    f_phys = [np.fft.ifftn(f_hat[i]).real for i in range(3)]
    rms = np.sqrt(np.mean(f_phys[0] ** 2 + f_phys[1] ** 2 + f_phys[2] ** 2))
    scale = amp / max(rms, 1e-12)
    return f_hat * scale


def main() -> int:
    parser = argparse.ArgumentParser(description="NS-BU-04 SGS frequency response")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out-dir", default=str(ROOT / "python" / "artifacts" / "paper"))
    args = parser.parse_args()

    N_hi = 32
    N_co = 16
    dt = 5e-4
    t_max = 0.6
    nu = 1e-3
    mu = 0.0
    alpha = 2.0
    dealias = True
    seed = args.seed
    k_init = 8
    record_u_hat_every = 10
    k_force = 2
    amp = 1.0e-2
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cfg = {
        "script": Path(__file__).name,
        "seed": seed,
        "N": N_hi,
        "N_co": N_co,
        "dt": dt,
        "t_max": t_max,
        "nu": nu,
        "mu": mu,
        "alpha": alpha,
        "k_force": k_force,
        "amp": amp,
        "k_init": k_init,
        "burn_in_frac": 0.4,
        "out_dir": str(out_dir),
    }
    print(f"FIGURE_CONFIG: {json.dumps(cfg, sort_keys=True)}")

    u0_hat = make_deterministic_u0_hat_3d(
        N_hi, L=2 * np.pi, seed=seed, k_init=k_init, dealias=dealias
    )
    forcing_hat_hi = make_lowk_divfree_forcing_hat_3d(
        N_hi, L=2 * np.pi, seed=seed, k_force=k_force, amp=amp, dealias=dealias
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

    burn_idx = int(0.4 * len(res.t_hat))
    t_use = res.t_hat[burn_idx:]
    u_hist = res.u_hat_hist[burn_idx:]

    kx_hi, ky_hi, kz_hi, k2_hi = kgrid_3d(N_hi)
    kx_co, ky_co, kz_co, k2_co = kgrid_3d(N_co)
    mask_hi = dealias_mask_3d(N_hi) if dealias else None
    mask_co = dealias_mask_3d(N_co) if dealias else None

    u_series = []
    c_series = []
    for u_hi in u_hist:
        u_co = restrict_hat_3d(u_hi, N_co)
        if mask_co is not None:
            u_co = u_co * mask_co

        nl_hi = nonlinear_term_hat_3d(u_hi, kx_hi, ky_hi, kz_hi, k2_hi, mask=mask_hi)
        nl_hi_co = restrict_hat_3d(nl_hi, N_co)
        nl_co = nonlinear_term_hat_3d(u_co, kx_co, ky_co, kz_co, k2_co, mask=mask_co)

        tau = nl_hi_co - nl_co
        c_hat = -tau

        u_series.append(u_co)
        c_series.append(c_hat)

    u_series = np.array(u_series)
    c_series = np.array(c_series)

    window = np.hanning(len(t_use))
    u_series = u_series * window[:, None, None, None, None]
    c_series = c_series * window[:, None, None, None, None]

    U = np.fft.fft(u_series, axis=0)
    C = np.fft.fft(c_series, axis=0)
    dt_sample = float(t_use[1] - t_use[0]) if len(t_use) > 1 else dt
    freqs = np.fft.fftfreq(len(t_use), d=dt_sample)
    omega = 2 * np.pi * freqs

    k_mag = np.sqrt(k2_co)
    shell_idx = dyadic_shell_index(k_mag, k0=1.0)
    j_max = int(shell_idx.max())
    shells = [j for j in range(1, min(3, j_max) + 1)]

    U_flat = U.reshape(U.shape[0], 3, -1)
    C_flat = C.reshape(C.shape[0], 3, -1)
    shell_flat = shell_idx.ravel()

    H_shell = {}
    denom_shell = {}
    passivity_rows = []
    denom_tol = 1e-12
    for s in shells:
        mask = shell_flat == s
        if not np.any(mask):
            continue
        numer = np.sum(C_flat[:, :, mask] * np.conj(U_flat[:, :, mask]), axis=(1, 2))
        denom = np.sum(np.abs(U_flat[:, :, mask]) ** 2, axis=(1, 2))
        H = numer / np.maximum(denom, denom_tol)
        H_shell[s] = H
        denom_shell[s] = denom

        valid = (denom > denom_tol) & (omega > 0)
        if np.any(valid):
            re_vals = np.real(H[valid])
            re_min = float(np.min(re_vals))
            re_med = float(np.median(re_vals))
            frac_nonneg = float(np.mean(re_vals >= -1e-3))
            passivity_rows.append((s, int(np.sum(valid)), re_min, re_med, frac_nonneg))
        else:
            passivity_rows.append((s, 0, float("nan"), float("nan"), float("nan")))

    print("NS-BU-04 SGS frequency response")
    print(f"N_hi={N_hi} N_co={N_co} dt={dt} t_max={t_max} nu={nu} mu={mu} alpha={alpha}")
    print(f"forcing: k_force={k_force} amp={amp}")
    print("sign convention: tau = NL_hi->co - NL_co; C = -tau (C appears on LHS)")
    print(f"burn_in_frac=0.4 n_samples={len(t_use)} j_max={j_max}")
    print("shell | n_freq | re_min | re_median | frac_nonneg")
    for s, n_freq, re_min, re_med, frac_nonneg in passivity_rows:
        print(f"{s:5d} | {n_freq:6d} | {re_min:7.4f} | {re_med:9.4f} | {frac_nonneg:11.4f}")
    print("NOTE: short-time forced run; H_s(omega) is diagnostic and resolution-limited.")

    # time-domain power proxy per shell
    u_flat_t = u_series.reshape(u_series.shape[0], 3, -1)
    c_flat_t = c_series.reshape(c_series.shape[0], 3, -1)
    print("shell | P_mean")
    for s in shells:
        mask = shell_flat == s
        if not np.any(mask):
            continue
        dot = np.real(np.sum(c_flat_t[:, :, mask] * np.conj(u_flat_t[:, :, mask]), axis=(1, 2)))
        P_mean = float(np.mean(dot))
        print(f"{s:5d} | {P_mean:8.6f}")

    plt.figure(figsize=(7, 4))
    for s in shells:
        H = H_shell.get(s)
        if H is None:
            continue
        pos = omega > 0
        plt.plot(omega[pos], np.real(H)[pos], label=f"Re H_s (s={s})")
    plt.axhline(0.0, color="k", linestyle="--", linewidth=0.8)
    plt.xlabel("omega")
    plt.ylabel("Re H_s(omega)")
    plt.title("NS-BU-04 SGS response (shellwise)")
    plt.legend(loc="best")
    plt.tight_layout()

    out_path = out_dir / "nsbu04_sgs_H.png"
    plt.savefig(out_path, dpi=150)
    print(f"saved {out_path}")

    # save data for NS-BU-05
    H_arr = np.stack([H_shell[s] for s in shells], axis=0)
    denom_arr = np.stack([denom_shell[s] for s in shells], axis=0)
    data_path = out_dir / "nsbu04_sgs_H_data.npz"
    meta_json = np.frombuffer(json.dumps(cfg, sort_keys=True).encode("utf-8"), dtype=np.uint8)
    savez_deterministic(
        data_path,
        meta_json=meta_json,
        omega=omega,
        H=H_arr,
        denom=denom_arr,
        shells=np.array(shells),
    )
    print(f"saved data: {data_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
