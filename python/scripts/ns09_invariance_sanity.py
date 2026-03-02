import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "python" / "src"
sys.path.insert(0, str(SRC))

from nswave.operators import kgrid_2d  # noqa: E402


def vorticity_to_velocity_hat(omega_hat, kx, ky, k2):
    psi_hat = np.zeros_like(omega_hat)
    mask = k2 != 0
    psi_hat[mask] = omega_hat[mask] / k2[mask]
    u_hat = np.empty((2,) + omega_hat.shape, dtype=complex)
    u_hat[0] = 1j * ky * psi_hat
    u_hat[1] = -1j * kx * psi_hat
    return u_hat


def rate_from_closure(u_hat, c_hat, L):
    u0 = np.fft.ifft2(u_hat[0]).real
    u1 = np.fft.ifft2(u_hat[1]).real
    c0 = np.fft.ifft2(c_hat[0]).real
    c1 = np.fft.ifft2(c_hat[1]).real
    return float(np.mean(u0 * c0 + u1 * c1) * L * L)


def main() -> None:
    L = 2 * np.pi
    N = 32
    M = 20
    nu0 = 1e-2
    gamma = 0.5

    kx, ky, k2 = kgrid_2d(N, L=L)

    rng = np.random.default_rng(0)
    rate_good = []
    rate_aniso_x = []
    rate_aniso_y = []
    rate_mean = []
    rate_anti = []

    for _ in range(M):
        omega = rng.standard_normal((N, N))
        omega_hat = np.fft.fft2(omega)
        u_hat = vorticity_to_velocity_hat(omega_hat, kx, ky, k2)

        c_good = -nu0 * k2 * u_hat
        c_good[:, 0, 0] = 0.0

        c_aniso_x = -nu0 * (kx * kx) * u_hat
        c_aniso_y = -nu0 * (ky * ky) * u_hat
        c_aniso_x[:, 0, 0] = 0.0
        c_aniso_y[:, 0, 0] = 0.0

        c_mean = np.zeros_like(u_hat)
        c_mean[:, 0, 0] = -gamma * u_hat[:, 0, 0]

        c_anti = nu0 * k2 * u_hat
        c_anti[:, 0, 0] = 0.0

        rate_good.append(rate_from_closure(u_hat, c_good, L))
        rate_aniso_x.append(rate_from_closure(u_hat, c_aniso_x, L))
        rate_aniso_y.append(rate_from_closure(u_hat, c_aniso_y, L))
        rate_mean.append(rate_from_closure(u_hat, c_mean, L))
        rate_anti.append(rate_from_closure(u_hat, c_anti, L))

    rate_good = np.array(rate_good)
    rate_aniso_x = np.array(rate_aniso_x)
    rate_aniso_y = np.array(rate_aniso_y)
    rate_mean = np.array(rate_mean)
    rate_anti = np.array(rate_anti)

    frac_good = float(np.mean(rate_good <= 1e-12))
    frac_anti = float(np.mean(rate_anti >= -1e-12))

    denom = np.maximum(np.maximum(np.abs(rate_aniso_x), np.abs(rate_aniso_y)), 1e-30)
    anis_ratio = np.abs(rate_aniso_x - rate_aniso_y) / denom
    anis_median = float(np.median(anis_ratio))

    u_const_x = np.ones((N, N))
    u_const_y = np.zeros((N, N))
    u_hat_const = np.stack([np.fft.fft2(u_const_x), np.fft.fft2(u_const_y)], axis=0)
    c_good_const = -nu0 * k2 * u_hat_const
    c_good_const[:, 0, 0] = 0.0
    c_mean_const = np.zeros_like(u_hat_const)
    c_mean_const[:, 0, 0] = -gamma * u_hat_const[:, 0, 0]

    rate_good_const = rate_from_closure(u_hat_const, c_good_const, L)
    rate_mean_const = rate_from_closure(u_hat_const, c_mean_const, L)

    print("NS-09 invariance sanity")
    print(f"N={N} M={M} nu0={nu0} gamma={gamma}")
    print(f"frac_good_passive={frac_good:.2f}")
    print(f"frac_anti_injects={frac_anti:.2f}")
    print(f"anis_ratio_median={anis_median:.3f}")
    print(f"rate_good_const={rate_good_const:.3e}")
    print(f"rate_mean_const={rate_mean_const:.3e}")

    print(f"BAD anisotropic: violates isotropy (anis_ratio median = {anis_median:.3f})")
    print(f"BAD mean damping: damps k=0 (rate_const = {rate_mean_const:.3e})")
    print(f"BAD anti-diffusion: injects energy (rate median = {np.median(rate_anti):.3e})")

    # Single-mode anisotropy check
    x = L * np.arange(N) / N
    y = L * np.arange(N) / N
    xx, yy = np.meshgrid(x, y, indexing="ij")
    omega_single = np.cos(6 * xx)
    omega_hat_single = np.fft.fft2(omega_single)
    u_hat_single = vorticity_to_velocity_hat(omega_hat_single, kx, ky, k2)
    c_aniso_x_single = -nu0 * (kx * kx) * u_hat_single
    c_aniso_y_single = -nu0 * (ky * ky) * u_hat_single
    c_aniso_x_single[:, 0, 0] = 0.0
    c_aniso_y_single[:, 0, 0] = 0.0
    rate_aniso_x_single = rate_from_closure(u_hat_single, c_aniso_x_single, L)
    rate_aniso_y_single = rate_from_closure(u_hat_single, c_aniso_y_single, L)
    anis_ratio_single = abs(rate_aniso_x_single - rate_aniso_y_single) / max(
        abs(rate_aniso_x_single), abs(rate_aniso_y_single), 1e-30
    )
    print(f"single_mode_anis_ratio={anis_ratio_single:.3e}")

    artifacts = ROOT / "python" / "artifacts"
    artifacts.mkdir(parents=True, exist_ok=True)
    out_path = artifacts / "ns09_invariance.png"

    labels = ["good", "aniso_x", "mean", "anti"]
    medians = [np.median(rate_good), np.median(rate_aniso_x), np.median(rate_mean), np.median(rate_anti)]

    plt.figure(figsize=(6, 4))
    plt.bar(labels, medians)
    plt.axhline(0.0, color="k", linewidth=0.5)
    plt.ylabel("median rate")
    plt.title("NS-09 invariance sanity")
    plt.tight_layout()
    plt.savefig(out_path)
    print(f"saved {out_path}")


if __name__ == "__main__":
    main()
