import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "python" / "src"
sys.path.insert(0, str(SRC))

from nswave.coarse_grain import restrict_hat_3d
from nswave.dyadic import dyadic_shell_index
from nswave.ns3d_spectral import make_deterministic_u0_hat_3d, simulate_ns3d, nonlinear_term_hat_3d
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
    N_hi = 32
    N_mid = 16
    N_co = 8
    dt = 5e-4
    t_max = 0.2
    nu = 1e-3
    mu = 0.0
    alpha = 2.0
    dealias = True
    seed = 0
    k_init = 3
    record_u_hat_every = 20
    burn_in_frac = 0.4
    k_force = 2
    amp = 1.0e-2

    ART.mkdir(parents=True, exist_ok=True)

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

    burn_idx = int(burn_in_frac * len(res.t_hat))
    u_hist = res.u_hat_hist[burn_idx:]

    kx_hi, ky_hi, kz_hi, k2_hi = kgrid_3d(N_hi)
    kx_mid, ky_mid, kz_mid, k2_mid = kgrid_3d(N_mid)
    kx_co, ky_co, kz_co, k2_co = kgrid_3d(N_co)

    mask_hi = dealias_mask_3d(N_hi) if dealias else None
    mask_mid = dealias_mask_3d(N_mid) if dealias else None
    mask_co = dealias_mask_3d(N_co) if dealias else None

    k_mag_co = np.sqrt(k2_co)
    shell_idx = dyadic_shell_index(k_mag_co, k0=1.0)
    j_max = int(shell_idx.max())

    denom = np.zeros(j_max + 1)
    num_direct = np.zeros(j_max + 1)
    num_two = np.zeros(j_max + 1)

    max_u_co_diff = 0.0

    for u_hi in u_hist:
        u_mid = restrict_hat_3d(u_hi, N_mid)
        u_co_one = restrict_hat_3d(u_hi, N_co)
        u_co_two = restrict_hat_3d(u_mid, N_co)

        if mask_mid is not None:
            u_mid *= mask_mid
        if mask_co is not None:
            u_co_one *= mask_co
            u_co_two *= mask_co

        max_u_co_diff = max(max_u_co_diff, float(np.max(np.abs(u_co_one - u_co_two))))

        nl_hi = nonlinear_term_hat_3d(u_hi, kx_hi, ky_hi, kz_hi, k2_hi, mask=mask_hi)
        nl_mid = nonlinear_term_hat_3d(u_mid, kx_mid, ky_mid, kz_mid, k2_mid, mask=mask_mid)
        nl_co = nonlinear_term_hat_3d(u_co_one, kx_co, ky_co, kz_co, k2_co, mask=mask_co)

        nl_hi_to_co = restrict_hat_3d(nl_hi, N_co)
        nl_mid_to_co = restrict_hat_3d(nl_mid, N_co)

        c_direct = nl_co - nl_hi_to_co
        c_two = nl_co - nl_mid_to_co

        shell_flat = shell_idx.ravel()
        u_flat = u_co_one.reshape(3, -1)
        c_direct_flat = c_direct.reshape(3, -1)
        c_two_flat = c_two.reshape(3, -1)

        for s in range(j_max + 1):
            mask = shell_flat == s
            if not np.any(mask):
                continue
            uu = np.sum(np.abs(u_flat[:, mask]) ** 2, axis=0)
            denom[s] += np.sum(uu)

            dot_direct = np.real(np.sum(c_direct_flat[:, mask] * np.conj(u_flat[:, mask]), axis=0))
            dot_two = np.real(np.sum(c_two_flat[:, mask] * np.conj(u_flat[:, mask]), axis=0))
            num_direct[s] += np.sum(dot_direct)
            num_two[s] += np.sum(dot_two)

    eps = 1e-12
    ell_direct = num_direct / np.maximum(denom, eps)
    ell_two = num_two / np.maximum(denom, eps)

    w = denom / max(np.sum(denom), eps)
    diff = ell_direct - ell_two
    rel_L2 = float(np.sqrt(np.sum(w * diff**2)) / max(np.sqrt(np.sum(w * ell_direct**2)), eps))
    max_abs = float(np.max(np.abs(diff)))
    max_rel = float(np.max(np.abs(diff) / (np.abs(ell_direct) + eps)))

    print("NS-BU-06 route dependence")
    print(f"N_hi={N_hi} N_mid={N_mid} N_co={N_co} dt={dt} t_max={t_max} nu={nu} mu={mu} alpha={alpha}")
    print(f"forcing: k_force={k_force} amp={amp} burn_in_frac={burn_in_frac}")
    print(f"j_max={j_max}")
    print(f"max_u_co_diff={max_u_co_diff:.6e}")
    print(f"rel_L2={rel_L2:.6e} max_abs={max_abs:.6e} max_rel={max_rel:.6e}")
    print(f"ell_direct={ell_direct}")
    print(f"ell_two={ell_two}")
    print(f"min_ell_direct={float(np.min(ell_direct)):.6e} min_ell_two={float(np.min(ell_two)):.6e}")

    fig, axes = plt.subplots(2, 1, figsize=(7, 7), sharex=True)
    shells = np.arange(j_max + 1)
    axes[0].plot(shells, ell_direct, marker="o", label="direct")
    axes[0].plot(shells, ell_two, marker="s", label="two-step")
    axes[0].set_ylabel("ell_s")
    axes[0].legend(loc="best")

    axes[1].plot(shells, np.abs(diff) / (np.abs(ell_direct) + eps), marker="o")
    axes[1].set_xlabel("shell s")
    axes[1].set_ylabel("|diff|/|ell_direct|")

    plt.tight_layout()
    out_path = ART / "nsbu06_route_dependence.png"
    plt.savefig(out_path, dpi=150)
    print(f"saved {out_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
