import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "python" / "src"
sys.path.insert(0, str(SRC))

from nswave.ns2d_spectral import make_deterministic_omega0_hat, simulate_ns2d  # noqa: E402
from nswave.operators import kgrid_2d  # noqa: E402


def low_k_error(
    omega_hat_base: np.ndarray,
    omega_hat_mu: np.ndarray,
    k2: np.ndarray,
    L: float,
    k_low: float,
) -> float:
    k_mag = np.sqrt(k2)
    mask = k_mag <= k_low

    omega_low_base = np.fft.ifft2(omega_hat_base * mask).real
    omega_low_mu = np.fft.ifft2(omega_hat_mu * mask).real

    norm_base = np.sqrt(np.mean(omega_low_base**2) * L * L)
    norm_err = np.sqrt(np.mean((omega_low_mu - omega_low_base) ** 2) * L * L)
    return float(norm_err / max(norm_base, 1e-12))


def main() -> None:
    N = 64
    L = 2 * np.pi
    dt = 5e-4
    t_max = 0.2
    nu = 1e-3
    alpha = 2.0
    mu_list = [0.0, 1e-2, 3e-3, 1e-3]
    dealias = True
    k_low = 6

    omega0_hat = make_deterministic_omega0_hat(N, dealias=dealias)
    kx, ky, k2 = kgrid_2d(N, L=L)

    base = simulate_ns2d(
        N=N,
        dt=dt,
        t_max=t_max,
        nu=nu,
        mu=0.0,
        alpha=alpha,
        omega0_hat=omega0_hat,
        forcing_hat_fn=None,
        L=L,
        dealias=dealias,
        record_every=10,
    )
    omega_hat_base = base.omega_hat_final

    errs = {}
    for mu in mu_list[1:]:
        res = simulate_ns2d(
            N=N,
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
        )
        errs[mu] = low_k_error(omega_hat_base, res.omega_hat_final, k2, L, k_low)

    print("NS-07 mu->0 convergence (low-k)")
    print(f"N={N} dt={dt} t_max={t_max} nu={nu} alpha={alpha} k_low={k_low}")
    print("mu | err_lowk")
    print("0.0 | 0.0")
    for mu in mu_list[1:]:
        print(f"{mu:.1e} | {errs[mu]:.6e}")

    mu_sorted = sorted(errs.keys(), reverse=True)
    monotone = all(errs[mu_sorted[i + 1]] <= errs[mu_sorted[i]] for i in range(len(mu_sorted) - 1))
    print(f"monotone_decreasing = {monotone}")

    artifacts = ROOT / "python" / "artifacts"
    artifacts.mkdir(parents=True, exist_ok=True)
    out_path = artifacts / "ns07_mu_convergence.png"

    mu_vals = np.array(mu_sorted)[::-1]
    err_vals = np.array([errs[mu] for mu in mu_vals])

    plt.figure(figsize=(5, 4))
    plt.loglog(mu_vals, err_vals, marker="o")
    plt.xlabel("mu")
    plt.ylabel("low-k relative error")
    plt.tight_layout()
    plt.savefig(out_path)
    print(f"saved {out_path}")


if __name__ == "__main__":
    main()
