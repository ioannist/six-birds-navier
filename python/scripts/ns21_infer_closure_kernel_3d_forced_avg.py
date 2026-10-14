import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "python" / "src"
sys.path.insert(0, str(SRC))

from nswave.coarse_grain import restrict_hat_3d, shell_index_from_k2  # noqa: E402
from nswave.ns3d_spectral import (  # noqa: E402
    make_deterministic_u0_hat_3d,
    nonlinear_term_hat_3d,
    simulate_ns3d,
)
from nswave.operators import dealias_mask_3d, kgrid_3d, project_div_free_3d  # noqa: E402


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
    low_mask = kmag2 <= k_force**2
    f_hat = f_hat * low_mask

    kx, ky, kz, k2 = kgrid_3d(N, L=L)
    f_hat = project_div_free_3d(f_hat, kx, ky, kz)
    if dealias:
        f_hat = f_hat * dealias_mask_3d(N)

    f_phys = [np.fft.ifftn(f_hat[i]).real for i in range(3)]
    rms = np.sqrt(np.mean(f_phys[0] ** 2 + f_phys[1] ** 2 + f_phys[2] ** 2))
    scale = amp / max(rms, 1e-12)
    return f_hat * scale


def accumulate_window(
    u_hat_hist: np.ndarray,
    t_hat: np.ndarray,
    t_start: float,
    t_end: float,
    *,
    N_hi: int,
    N_co: int,
    L: float,
    nu: float,
    dealias: bool,
    f_hat_hi: np.ndarray,
    f_hat_co: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    kx_hi, ky_hi, kz_hi, k2_hi = kgrid_3d(N_hi, L=L)
    mask_hi = dealias_mask_3d(N_hi) if dealias else None

    kx, ky, kz, k2 = kgrid_3d(N_co, L=L)
    mask_co = dealias_mask_3d(N_co) if dealias else None

    shell = shell_index_from_k2(k2)
    s_max = int(shell.max())
    num = np.zeros(s_max + 1, dtype=float)
    den = np.zeros(s_max + 1, dtype=float)

    for i in range(len(t_hat) - 1):
        t = t_hat[i]
        t_next = t_hat[i + 1]
        if t < t_start or t_next > t_end:
            continue
        dt_s = t_next - t
        if dt_s <= 0:
            continue

        u_hi = u_hat_hist[i]
        u_hi_next = u_hat_hist[i + 1]
        if mask_hi is not None:
            u_hi = u_hi * mask_hi
            u_hi_next = u_hi_next * mask_hi

        u_co = restrict_hat_3d(u_hi, N_co)
        u_co = project_div_free_3d(u_co, kx, ky, kz)
        u_co_next = restrict_hat_3d(u_hi_next, N_co)
        u_co_next = project_div_free_3d(u_co_next, kx, ky, kz)
        if mask_co is not None:
            u_co = u_co * mask_co
            u_co_next = u_co_next * mask_co

        w_t = (u_co_next - u_co) / dt_s

        lin_co = -(nu * k2)
        nl_co = nonlinear_term_hat_3d(
            u_co, kx, ky, kz, k2, mask=mask_co, t=t, forcing_hat_fn=None
        )
        rhs_no = lin_co * u_co + nl_co + f_hat_co

        tau = w_t - rhs_no
        if mask_co is not None:
            tau = tau * mask_co

        dot = np.real(
            tau[0] * np.conj(u_co[0])
            + tau[1] * np.conj(u_co[1])
            + tau[2] * np.conj(u_co[2])
        )
        uu = np.abs(u_co[0]) ** 2 + np.abs(u_co[1]) ** 2 + np.abs(u_co[2]) ** 2

        num += np.bincount(shell.ravel(), weights=-dot.ravel(), minlength=s_max + 1)
        den += np.bincount(shell.ravel(), weights=uu.ravel(), minlength=s_max + 1)

    return num, den


def fit_low_band(ell: np.ndarray, den: np.ndarray, s_lo: int, s_hi: int) -> tuple[float, float]:
    shells = np.arange(len(ell))
    mask = (shells >= s_lo) & (shells <= s_hi) & (den > 0)
    if not np.any(mask):
        return float("nan"), float("nan")
    s = shells[mask]
    ell_s = ell[mask]
    weights = den[mask]
    nu_eff = float(np.sum(weights * (s**2) * ell_s) / np.sum(weights * (s**4)))
    fit = nu_eff * s**2
    rel_err = float(
        np.sqrt(np.average((ell_s - fit) ** 2, weights=weights))
        / max(np.average(np.abs(ell_s), weights=weights), 1e-12)
    )
    return nu_eff, rel_err


def high_band_slope(ell: np.ndarray, den: np.ndarray, s_lo: int, s_hi: int) -> tuple[float, int]:
    shells = np.arange(len(ell))
    den_thresh = 1e-6 * float(np.max(den)) if np.max(den) > 0 else 0.0
    ell_abs = np.abs(ell)
    valid = (shells >= s_lo) & (shells <= s_hi) & (den >= den_thresh) & (ell_abs > 1e-14)
    slopes = []
    for s in range(s_lo, s_hi):
        if valid[s] and valid[s + 1]:
            slope = (np.log(ell_abs[s + 1]) - np.log(ell_abs[s])) / (np.log(s + 1) - np.log(s))
            slopes.append(float(slope))
    if slopes:
        return float(np.median(slopes)), len(slopes)
    return float("nan"), 0


def main() -> int:
    L = 2 * np.pi
    N_hi = 32
    N_co = 16
    dt = 5e-4
    t_max = 0.30
    nu = 1e-3
    mu = 0.0
    alpha = 2.0
    dealias = True
    seed_u0 = 0
    seed_force = 0
    k_init = 10
    k_force = 2
    amp = 1e-2
    record_u_hat_every = 10
    burn_in = 0.15
    n_windows = 3
    transient_window = (0.0, 0.02)

    ART.mkdir(parents=True, exist_ok=True)

    u0_hat = make_deterministic_u0_hat_3d(
        N_hi, L=L, seed=seed_u0, k_init=k_init, dealias=dealias
    )
    f_hat_hi = make_lowk_divfree_forcing_hat_3d(
        N_hi, L=L, seed=seed_force, k_force=k_force, amp=amp, dealias=dealias
    )
    f_hat_co = restrict_hat_3d(f_hat_hi, N_co)
    kx, ky, kz, _ = kgrid_3d(N_co, L=L)
    f_hat_co = project_div_free_3d(f_hat_co, kx, ky, kz)

    def forcing_hat_fn(_t: float, _kx: np.ndarray, _ky: np.ndarray, _kz: np.ndarray) -> np.ndarray:
        return f_hat_hi

    res_hi = simulate_ns3d(
        N=N_hi,
        dt=dt,
        t_max=t_max,
        nu=nu,
        mu=mu,
        alpha=alpha,
        u0_hat=u0_hat,
        forcing_hat_fn=forcing_hat_fn,
        L=L,
        dealias=dealias,
        record_every=10,
        record_u_hat_every=record_u_hat_every,
    )

    if res_hi.u_hat_hist is None or res_hi.t_hat is None:
        raise RuntimeError("u_hat_hist or t_hat missing in res_hi")

    t_hat = res_hi.t_hat
    u_hat_hist = res_hi.u_hat_hist

    t_start = burn_in
    t_end = t_max
    window_len = (t_end - t_start) / n_windows
    windows = [(t_start + i * window_len, t_start + (i + 1) * window_len) for i in range(n_windows)]

    num_list = []
    den_list = []
    for w_start, w_end in windows:
        num, den = accumulate_window(
            u_hat_hist,
            t_hat,
            w_start,
            w_end,
            N_hi=N_hi,
            N_co=N_co,
            L=L,
            nu=nu,
            dealias=dealias,
            f_hat_hi=f_hat_hi,
            f_hat_co=f_hat_co,
        )
        num_list.append(num)
        den_list.append(den)

    num_arr = np.array(num_list)
    den_arr = np.array(den_list)
    ell_list = num_arr / np.maximum(den_arr, 1e-30)

    ell_mean = np.mean(ell_list, axis=0)
    ell_lo = np.min(ell_list, axis=0)
    ell_hi = np.max(ell_list, axis=0)
    den_mean = np.mean(den_arr, axis=0)

    num_trans, den_trans = accumulate_window(
        u_hat_hist,
        t_hat,
        transient_window[0],
        transient_window[1],
        N_hi=N_hi,
        N_co=N_co,
        L=L,
        nu=nu,
        dealias=dealias,
        f_hat_hi=f_hat_hi,
        f_hat_co=f_hat_co,
    )
    ell_trans = num_trans / np.maximum(den_trans, 1e-30)

    ell_mean_abs = np.abs(ell_mean)
    ell_trans_abs = np.abs(ell_trans)

    nu_eff, rel_err_steady = fit_low_band(ell_mean_abs, den_mean, 2, 4)
    _, rel_err_trans = fit_low_band(ell_trans_abs, den_trans, 2, 4)

    shells = np.arange(len(ell_mean_abs))
    den_pos = den_mean > 0
    s_nonzero_max = int(shells[den_pos].max()) if np.any(den_pos) else 0
    r = np.zeros_like(ell_mean_abs)
    mask_r = den_pos & (shells >= 1)
    r[mask_r] = ell_mean_abs[mask_r] / np.maximum(nu_eff * (shells[mask_r] ** 2), 1e-30)
    low_mask_r = mask_r & (shells <= 4)
    high_start = int(np.ceil(0.7 * s_nonzero_max)) if s_nonzero_max > 0 else 0
    high_mask_r = mask_r & (shells >= high_start)
    steep_hi = float(np.median(r[high_mask_r])) if np.any(high_mask_r) else float("nan")

    diss_shell = np.sum(num_arr, axis=0)
    diss_total = float(np.sum(diss_shell[1:]))
    frac_high = (
        float(np.sum(diss_shell[high_start:])) / max(diss_total, 1e-30)
        if high_start > 0
        else float("nan")
    )

    high_band_p, n_pairs = high_band_slope(ell_mean_abs, den_mean, 6, 10)

    dissip_total = float(np.sum(np.sum(num_arr[:, 1:], axis=1)))
    closure_dissipative = dissip_total >= -1e-8

    mask_energy = res_hi.t >= burn_in
    if np.any(mask_energy):
        e_window = res_hi.energy[mask_energy]
        energy_drift = float((np.max(e_window) - np.min(e_window)) / max(np.mean(e_window), 1e-12))
    else:
        energy_drift = float("nan")

    print("NS-21 closure kernel inference (3D forced + averaged)")
    print(f"N_hi={N_hi} N_co={N_co} dt={dt} t_max={t_max} nu={nu} mu={mu}")
    print(f"forcing: k_force={k_force} amp={amp}")
    print(f"windows: burn_in={burn_in} n_windows={n_windows}")
    print(
        "nu_eff={:.6e} low_band_rel_err_steady={:.6e} low_band_rel_err_transient={:.6e}".format(
            nu_eff, rel_err_steady, rel_err_trans
        )
    )
    print(f"steep_hi={steep_hi:.3e}")
    print(f"frac_high={frac_high:.3e} high_start={high_start} s_nonzero_max={s_nonzero_max}")
    if n_pairs >= 6 and np.isfinite(high_band_p):
        print(f"high_band_p={high_band_p:.3f} n_pairs={n_pairs}")
        if high_band_p <= 2.0 or high_band_p > 10.0:
            print("WARNING: high_band_p outside expected range; fit likely noisy")
    else:
        print(f"high_band_p=SKIP n_pairs={n_pairs} (insufficient)")
    print(f"dissip_total={dissip_total:.6e} closure_dissipative={closure_dissipative}")
    print(f"energy_drift_window={energy_drift:.6e}")

    shells = np.arange(len(ell_mean))
    mask_plot = (shells >= 1) & (den_mean > 0) & (ell_mean > 0)
    s_plot = shells[mask_plot]
    ell_plot = ell_mean[mask_plot]

    plt.figure(figsize=(6, 4))
    plt.loglog(s_plot, ell_plot, "o-", label="mean")
    if np.any(mask_plot):
        lo = np.maximum(ell_lo[mask_plot], 1e-20)
        hi = np.maximum(ell_hi[mask_plot], 1e-20)
        plt.fill_between(s_plot, lo, hi, alpha=0.2, label="window band")
    if np.isfinite(nu_eff):
        s_low = s_plot[(s_plot >= 1) & (s_plot <= 4)]
        plt.loglog(s_low, nu_eff * s_low**2, "--", label="low-k fit")
    if np.isfinite(high_band_p) and n_pairs >= 6:
        s_hi = s_plot[(s_plot >= 6) & (s_plot <= 10)]
        if s_hi.size > 0:
            c_fit = np.exp(np.median(np.log(ell_plot[(s_plot >= 6) & (s_plot <= 10)]) - high_band_p * np.log(s_hi)))
            plt.loglog(s_hi, c_fit * s_hi**high_band_p, ":", label="high-k fit")
    plt.xlabel("shell s")
    plt.ylabel("ell(s)")
    plt.title(f"NS-21 closure kernel (N_hi={N_hi}→N_co={N_co})")
    plt.legend()
    plt.tight_layout()

    out_path = ART / "ns21_closure_kernel_3d_forced_avg.png"
    plt.savefig(out_path, dpi=150)
    print(f"saved {out_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
