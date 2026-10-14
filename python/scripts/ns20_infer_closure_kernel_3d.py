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


def main() -> int:
    L = 2 * np.pi
    N_hi = 32
    N_co = 16
    dt = 5e-4
    t_max = 0.05
    nu = 1e-3
    mu = 0.0
    alpha = 2.0
    dealias = True
    record_u_hat_every = 10
    t_min = 0.015
    k_init = 7

    ART.mkdir(parents=True, exist_ok=True)

    u0_hat = make_deterministic_u0_hat_3d(
        N_hi, L=L, seed=0, k_init=k_init, dealias=dealias
    )

    res_hi = simulate_ns3d(
        N=N_hi,
        dt=dt,
        t_max=t_max,
        nu=nu,
        mu=mu,
        alpha=alpha,
        u0_hat=u0_hat,
        forcing_hat_fn=None,
        L=L,
        dealias=dealias,
        record_every=10,
        record_u_hat_every=record_u_hat_every,
    )

    if res_hi.u_hat_hist is None or res_hi.t_hat is None:
        raise RuntimeError("u_hat_hist or t_hat missing in res_hi")

    kx_hi, ky_hi, kz_hi, k2_hi = kgrid_3d(N_hi, L=L)
    mask_hi = dealias_mask_3d(N_hi) if dealias else None

    kx, ky, kz, k2 = kgrid_3d(N_co, L=L)
    mask_co = dealias_mask_3d(N_co) if dealias else None

    shell = shell_index_from_k2(k2)
    s_max = int(shell.max())
    num = np.zeros(s_max + 1, dtype=float)
    den = np.zeros(s_max + 1, dtype=float)

    for t, u_hi in zip(res_hi.t_hat, res_hi.u_hat_hist):
        if t < t_min:
            continue
        if mask_hi is not None:
            u_hi = u_hi * mask_hi

        u_co = restrict_hat_3d(u_hi, N_co)
        u_co = project_div_free_3d(u_co, kx, ky, kz)
        if mask_co is not None:
            u_co = u_co * mask_co

        lin_hi = -(nu * k2_hi)
        nl_hi = nonlinear_term_hat_3d(
            u_hi, kx_hi, ky_hi, kz_hi, k2_hi, mask=mask_hi, t=t, forcing_hat_fn=None
        )
        rhs_hi = lin_hi * u_hi + nl_hi
        rhs_hi_co = restrict_hat_3d(rhs_hi, N_co)
        rhs_hi_co = project_div_free_3d(rhs_hi_co, kx, ky, kz)

        lin_co = -(nu * k2)
        nl_co = nonlinear_term_hat_3d(
            u_co, kx, ky, kz, k2, mask=mask_co, t=t, forcing_hat_fn=None
        )
        rhs_no = lin_co * u_co + nl_co

        tau = rhs_hi_co - rhs_no
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

    eps = 1e-30
    ell = num / np.maximum(den, eps)

    dissip_total = float(np.sum(num[1:]))
    closure_dissipative = dissip_total >= -1e-8

    shells = np.arange(len(ell))
    low_mask = (shells >= 1) & (shells <= 4) & (den > 0)
    if np.any(low_mask):
        s_low = shells[low_mask]
        ell_low = ell[low_mask]
        nu_eff = float(np.sum(s_low**2 * ell_low) / np.sum(s_low**4))
        fit_low = nu_eff * s_low**2
        rel_err = float(
            np.sqrt(np.mean((ell_low - fit_low) ** 2)) / max(np.mean(np.abs(ell_low)), 1e-12)
        )
    else:
        nu_eff = float("nan")
        rel_err = float("nan")

    den_pos = den > 0
    s_nonzero_max = int(shells[den_pos].max()) if np.any(den_pos) else 0
    ell_abs = np.abs(ell)
    r = np.zeros_like(ell_abs)
    mask_r = den_pos & (shells >= 1)
    r[mask_r] = ell_abs[mask_r] / np.maximum(nu_eff * (shells[mask_r] ** 2), 1e-30)
    low_mask_r = mask_r & (shells <= 4)
    high_start = int(np.ceil(0.7 * s_nonzero_max)) if s_nonzero_max > 0 else 0
    high_mask_r = mask_r & (shells >= high_start)
    steep_hi = float(np.median(r[high_mask_r])) if np.any(high_mask_r) else float("nan")

    diss_shell = num.copy()
    diss_total = float(np.sum(diss_shell[1:]))
    frac_high = (
        float(np.sum(diss_shell[high_start:])) / max(diss_total, 1e-30)
        if high_start > 0
        else float("nan")
    )

    s_hi_min = max(4, int(s_max * 0.5))
    s_hi_max = min(s_max, s_hi_min + 4)
    den_thresh = 1e-6 * float(np.max(den)) if np.max(den) > 0 else 0.0
    high_mask = (shells >= s_hi_min) & (shells <= s_hi_max) & (den >= den_thresh) & (ell_abs > 1e-14)
    slope_list = []
    for s in range(s_hi_min, s_hi_max):
        if high_mask[s] and high_mask[s + 1]:
            slope = (np.log(ell_abs[s + 1]) - np.log(ell_abs[s])) / (np.log(s + 1) - np.log(s))
            slope_list.append(float(slope))
    if slope_list:
        high_band_p = float(np.median(slope_list))
        c_fit = float(np.exp(np.median(np.log(ell_abs[high_mask]) - high_band_p * np.log(shells[high_mask]))))
        n_pairs_used = len(slope_list)
    else:
        high_band_p = float("nan")
        c_fit = float("nan")
        n_pairs_used = 0

    print("NS-20 closure kernel inference (3D)")
    print(f"N_hi={N_hi} N_co={N_co} dt={dt} t_max={t_max} nu={nu} mu={mu}")
    print(f"nu_eff={nu_eff:.6e} low_band_rel_err={rel_err:.6e}")
    print(f"steep_hi={steep_hi:.3e}")
    print(f"frac_high={frac_high:.3e} high_start={high_start} s_nonzero_max={s_nonzero_max}")

    if n_pairs_used >= 6 and np.isfinite(high_band_p):
        print(f"high_band_p={high_band_p:.3f} n_pairs={n_pairs_used} band=[{s_hi_min},{s_hi_max}]")
        if high_band_p <= 2.0:
            print("WARNING: high_band_p <= 2.0; fit may be noisy")
        if high_band_p > 10.0:
            print("WARNING: high_band_p > 10; fit likely noisy")
    else:
        print(f"high_band_p=SKIP n_pairs={n_pairs_used} band=[{s_hi_min},{s_hi_max}] (insufficient)")
    print(f"dissip_total={dissip_total:.6e} closure_dissipative={closure_dissipative}")

    mask_plot = (shells >= 1) & (den > 0)
    s_plot = shells[mask_plot]
    ell_plot = ell[mask_plot]

    plt.figure(figsize=(6, 4))
    plt.loglog(s_plot, ell_plot, "o-", label="inferred")
    if np.any(low_mask):
        plt.loglog(s_low, nu_eff * s_low**2, "--", label="low-k fit")
    if np.isfinite(high_band_p) and n_pairs_used >= 6:
        s_fit = s_plot[s_plot >= s_hi_min]
        if s_fit.size > 0:
            plt.loglog(s_fit, c_fit * s_fit**high_band_p, ":", label="high-k fit")
    plt.xlabel("shell s")
    plt.ylabel("ell(s)")
    plt.title(f"NS-20 closure kernel (N_hi={N_hi}→N_co={N_co})")
    plt.legend()
    plt.tight_layout()

    out_path = ART / "ns20_closure_kernel_3d.png"
    plt.savefig(out_path, dpi=150)
    print(f"saved {out_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
