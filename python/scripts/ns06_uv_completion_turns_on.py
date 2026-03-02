import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "python" / "src"
sys.path.insert(0, str(SRC))

from nswave.ns2d_spectral import (  # noqa: E402
    dealias_mask_2d,
    make_deterministic_omega0_hat,
    simulate_ns2d,
)
from nswave.operators import kgrid_2d  # noqa: E402


def initial_omega_hat(N: int, mask: np.ndarray | None) -> np.ndarray:
    if mask is None:
        return make_deterministic_omega0_hat(N, dealias=False)
    return make_deterministic_omega0_hat(N, dealias=True)


def energy_spectrum_shells(omega_hat: np.ndarray, k2: np.ndarray) -> np.ndarray:
    k_mag = np.sqrt(k2)
    shell = np.floor(k_mag + 0.5).astype(int)

    E_mode = np.zeros_like(k2, dtype=float)
    mask = k2 > 0
    E_mode[mask] = 0.5 * (np.abs(omega_hat[mask]) ** 2) / k2[mask]

    kmax = int(shell.max())
    E_shell = np.bincount(shell.ravel(), weights=E_mode.ravel(), minlength=kmax + 1)
    return E_shell


def main() -> None:
    N = 64
    L = 2 * np.pi
    dt = 5e-4
    t_max = 0.4
    nu = 1e-3
    alpha = 2.0
    mu_baseline = 0.0
    mu_completed = 1e-4
    record_omega_hat_every = int(round(0.05 / dt))
    dealias = True

    kx, ky, k2 = kgrid_2d(N, L=L)
    mask = dealias_mask_2d(N) if dealias else None
    omega0_hat = initial_omega_hat(N, mask)

    base = simulate_ns2d(
        N=N,
        dt=dt,
        t_max=t_max,
        nu=nu,
        mu=mu_baseline,
        alpha=alpha,
        omega0_hat=omega0_hat,
        forcing_hat_fn=None,
        L=L,
        dealias=dealias,
        record_every=10,
        record_omega_hat_every=record_omega_hat_every,
    )

    comp = simulate_ns2d(
        N=N,
        dt=dt,
        t_max=t_max,
        nu=nu,
        mu=mu_completed,
        alpha=alpha,
        omega0_hat=omega0_hat,
        forcing_hat_fn=None,
        L=L,
        dealias=dealias,
        record_every=10,
        record_omega_hat_every=record_omega_hat_every,
    )

    if base.omega_hat_hist is None or comp.omega_hat_hist is None:
        raise RuntimeError("omega_hat history was not recorded")

    t_targets = [0.1, 0.2, 0.4]
    k_low = 6
    k0 = 12

    print("NS-06 UV completion stress test")
    print(
        f"N={N} dt={dt} t_max={t_max} nu={nu} alpha={alpha} "
        f"mu_base={mu_baseline} mu_comp={mu_completed} k0={k0}"
    )

    for t_target in t_targets:
        idx = int(np.argmin(np.abs(base.t_omega - t_target)))
        t = float(base.t_omega[idx])

        E_base = energy_spectrum_shells(base.omega_hat_hist[idx], k2)
        E_comp = energy_spectrum_shells(comp.omega_hat_hist[idx], k2)

        low_base = float(np.sum(E_base[: k_low + 1]))
        low_comp = float(np.sum(E_comp[: k_low + 1]))
        high_base = float(np.sum(E_base[k0:]))
        high_comp = float(np.sum(E_comp[k0:]))

        ratio_low = low_comp / max(low_base, 1e-30)
        ratio_high = high_comp / max(high_base, 1e-30)

        print(
            f"t={t:.3f} ratio_low={ratio_low:.3f} ratio_high={ratio_high:.3f} "
            f"tail_base={high_base:.3e} tail_comp={high_comp:.3e}"
        )

    # Plot final time spectrum
    E_base = energy_spectrum_shells(base.omega_hat_hist[-1], k2)
    E_comp = energy_spectrum_shells(comp.omega_hat_hist[-1], k2)
    k = np.arange(len(E_base))

    artifacts = ROOT / "python" / "artifacts"
    artifacts.mkdir(parents=True, exist_ok=True)
    out_path = artifacts / "ns06_uv_completion.png"

    plt.figure(figsize=(6, 4))
    plt.semilogy(k[1:], E_base[1:], label="baseline (mu=0)")
    plt.semilogy(k[1:], E_comp[1:], label="completed (mu>0)")
    plt.xlabel("k")
    plt.ylabel("E(k)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path)
    print(f"saved {out_path}")


if __name__ == "__main__":
    main()
