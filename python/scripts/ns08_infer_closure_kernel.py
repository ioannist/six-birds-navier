import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "python" / "src"
sys.path.insert(0, str(SRC))

from nswave.coarse_grain import restrict_hat_2d, shell_index_from_k2  # noqa: E402
from nswave.ns2d_spectral import (
    dealias_mask_2d,
    make_deterministic_omega0_hat,
    nonlinear_term_hat,
    simulate_ns2d,
)  # noqa: E402
from nswave.operators import kgrid_2d  # noqa: E402


def main() -> None:
    L = 2 * np.pi
    N_hi = 128
    N_co = 64
    dt = 2e-4
    t_max = 0.20
    nu = 1e-3
    mu = 0.0
    alpha = 2.0
    dealias = True
    record_omega_hat_every = 20

    k_init = 40
    omega0_hat = make_deterministic_omega0_hat(N_hi, k_init=k_init, dealias=dealias)

    res_hi = simulate_ns2d(
        N=N_hi,
        dt=dt,
        t_max=t_max,
        nu=nu,
        mu=mu,
        alpha=alpha,
        omega0_hat=omega0_hat,
        forcing_hat_fn=None,
        L=L,
        dealias=dealias,
        record_every=10,
        record_omega_hat_every=record_omega_hat_every,
    )

    if res_hi.omega_hat_hist is None or res_hi.t_omega is None:
        raise RuntimeError("omega_hat history not recorded")

    t_hist = res_hi.t_omega
    omega_hi_hist = res_hi.omega_hat_hist

    kx_hi_c, ky_hi_c, k2_hi = kgrid_2d(N_hi, L=L)
    mask_hi = dealias_mask_2d(N_hi) if dealias else None

    kx, ky, k2 = kgrid_2d(N_co, L=L)
    mask_co = dealias_mask_2d(N_co) if dealias else None
    shell = shell_index_from_k2(k2)
    s_max = int(shell.max())

    num = np.zeros(s_max + 1)
    den = np.zeros(s_max + 1)
    tau_ratio_sum = 0.0
    tau_ratio_count = 0

    linop = -(nu * k2)
    linop_hi = -(nu * k2_hi)

    t_min = 0.02
    for i in range(len(t_hist)):
        if t_hist[i] < t_min:
            continue
        omega_hi = omega_hi_hist[i]
        if mask_hi is not None:
            omega_hi = omega_hi * mask_hi

        w_i = restrict_hat_2d(omega_hi, N_co)
        if mask_co is not None:
            w_i = w_i * mask_co

        nonlinear_hat_hi = nonlinear_term_hat(
            omega_hi, kx_hi_c, ky_hi_c, k2_hi, mask=mask_hi, t=t_hist[i]
        )
        rhs_hi = linop_hi * omega_hi + nonlinear_hat_hi
        rhs_hi_co = restrict_hat_2d(rhs_hi, N_co)

        nonlinear_hat = nonlinear_term_hat(w_i, kx, ky, k2, mask=mask_co, t=t_hist[i])
        rhs_no_closure = linop * w_i + nonlinear_hat
        tau = rhs_hi_co - rhs_no_closure
        if mask_co is not None:
            tau = tau * mask_co

        max_tau = np.max(np.abs(tau))
        max_w = np.max(np.abs(w_i))
        tau_ratio_sum += max_tau / max(max_w, 1e-30)
        tau_ratio_count += 1

        num += -np.bincount(
            shell.ravel(),
            weights=np.real(tau * np.conj(w_i)).ravel(),
            minlength=s_max + 1,
        )
        den += np.bincount(
            shell.ravel(), weights=np.abs(w_i).ravel() ** 2, minlength=s_max + 1
        )

    eps = 1e-30
    ell = num / np.maximum(den, eps)

    dissip_total = float(np.sum(num[1:]))
    closure_dissipative = dissip_total >= -1e-10

    low_band = (
        (np.arange(len(ell)) >= 1)
        & (np.arange(len(ell)) <= 6)
        & (den > 0)
        & (ell > 0)
    )
    if np.any(low_band):
        s_low = np.arange(len(ell))[low_band].astype(float)
        ell_low = ell[low_band]
        nu_eff = float(np.sum((s_low**2) * ell_low) / np.sum(s_low**4))
        fit = nu_eff * (s_low**2)
        rel_err = float(np.sqrt(np.mean((ell_low - fit) ** 2)) / max(np.mean(np.abs(ell_low)), 1e-12))
    else:
        nu_eff = float("nan")
        rel_err = float("nan")

    high_band = (np.arange(len(ell)) >= 10) & (np.arange(len(ell)) <= 18) & (ell > 1e-12)
    if np.count_nonzero(high_band) >= 2:
        s_hi = np.arange(len(ell))[high_band].astype(float)
        ell_hi = ell[high_band]
        coeffs = np.polyfit(np.log(s_hi), np.log(ell_hi), deg=1)
        p = float(coeffs[0])
        c = float(np.exp(coeffs[1]))
    else:
        p = None
        c = None

    ell_low_band = (np.arange(len(ell)) >= 2) & (np.arange(len(ell)) <= 6) & (den > 0)
    ell_high_band = (np.arange(len(ell)) >= 12) & (np.arange(len(ell)) <= 18) & (den > 0)
    if np.any(ell_low_band):
        s_low_m = np.arange(len(ell))[ell_low_band].astype(float)
        ell_low_median = float(np.median(ell[ell_low_band] / (s_low_m**2)))
    else:
        ell_low_median = float("nan")
    if np.any(ell_high_band):
        ell_high_median = float(np.median(ell[ell_high_band]))
    else:
        ell_high_median = float("nan")

    avg_tau_ratio = tau_ratio_sum / max(tau_ratio_count, 1)

    print("NS-08 closure kernel inference")
    print(f"N_hi={N_hi} N_co={N_co} dt={dt} t_max={t_max} nu={nu} mu={mu}")
    print(f"nu_eff={nu_eff:.6e} low_band_rel_err={rel_err:.3e}")
    if p is None:
        print("high_band_fit: insufficient points")
    else:
        print(f"high_band_p={p:.3f}")
        if p > 10:
            print("high_band_fit: p is very large; fit may be noisy")
    print(f"dissip_total={dissip_total:.6e} closure_dissipative={closure_dissipative}")
    print(f"ell_low_median={ell_low_median:.6e} ell_high_median={ell_high_median:.6e}")
    print(f"avg_max_abs_tau_over_w={avg_tau_ratio:.6e}")
    if ell_low_median < 1e-8 and dissip_total < 1e-10:
        print("WARNING: closure signal is tiny; increase k_init or t_max")

    artifacts = ROOT / "python" / "artifacts"
    artifacts.mkdir(parents=True, exist_ok=True)
    out_path = artifacts / "ns08_closure_kernel.png"

    s = np.arange(len(ell))
    valid = (den > 0) & (s >= 1)

    plt.figure(figsize=(6, 4))
    plt.loglog(s[valid], ell[valid], marker="o", linestyle="none", label="ell(s)")
    if np.any(low_band):
        s_line = s_low
        plt.loglog(s_line, nu_eff * s_line**2, label="low-k fit")
    if p is not None and c is not None:
        s_line = s[high_band]
        plt.loglog(s_line, c * s_line**p, label="high-k fit")
    plt.xlabel("shell s")
    plt.ylabel("ell(s)")
    plt.title(f"NS-08 closure kernel (N_hi={N_hi} -> N_co={N_co}, nu={nu})")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path)
    print(f"saved {out_path}")


if __name__ == "__main__":
    main()
