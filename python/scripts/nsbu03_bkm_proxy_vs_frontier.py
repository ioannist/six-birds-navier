import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "python" / "src"
sys.path.insert(0, str(SRC))

from nswave.dyadic import (
    crossing_times,
    frontier_scale,
    monotone_envelope,
    percentile_frontier,
    shell_energies,
)
from nswave.ns3d_spectral import make_deterministic_u0_hat_3d, simulate_ns3d
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


def vorticity_hat(u_hat: np.ndarray, kx: np.ndarray, ky: np.ndarray, kz: np.ndarray) -> np.ndarray:
    omega_hat = np.zeros_like(u_hat)
    omega_hat[0] = 1j * (ky * u_hat[2] - kz * u_hat[1])
    omega_hat[1] = 1j * (kz * u_hat[0] - kx * u_hat[2])
    omega_hat[2] = 1j * (kx * u_hat[1] - ky * u_hat[0])
    return omega_hat


def main() -> int:
    N = 32
    L = 2 * np.pi
    dt = 5e-4
    t_max = 0.3
    nu = 1e-3
    mu = 0.0
    alpha = 2.0
    seed = 0
    k_init = 3
    record_u_hat_every = 20
    k_force = 2
    amp = 1.0e-2

    ART.mkdir(parents=True, exist_ok=True)

    u0_hat = make_deterministic_u0_hat_3d(N, L=L, seed=seed, k_init=k_init, dealias=True)
    forcing_hat = make_lowk_divfree_forcing_hat_3d(
        N,
        L=L,
        seed=seed,
        k_force=k_force,
        amp=amp,
        dealias=True,
    )

    def forcing_hat_fn(t: float, kx: np.ndarray, ky: np.ndarray, kz: np.ndarray) -> np.ndarray:
        return forcing_hat

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
        record_every=10,
        record_u_hat_every=record_u_hat_every,
    )

    if res.u_hat_hist is None or res.t_hat is None:
        raise RuntimeError("u_hat_hist not recorded")

    kx, ky, kz, k2 = kgrid_3d(N, L=L)
    k_mag = np.sqrt(k2)

    omega_inf = []
    total_energy = []
    E_shell_all = []
    for u_hat in res.u_hat_hist:
        omega_hat = vorticity_hat(u_hat, kx, ky, kz)
        omega = [np.fft.ifftn(omega_hat[i]).real for i in range(3)]
        omega_mag = np.sqrt(omega[0] ** 2 + omega[1] ** 2 + omega[2] ** 2)
        omega_inf.append(float(np.max(omega_mag)))

        e_shell = shell_energies(u_hat, k_mag, L=L)
        E_shell_all.append(e_shell)
        total_energy.append(float(np.sum(e_shell)))

    omega_inf = np.array(omega_inf)
    total_energy = np.array(total_energy)
    t = np.array(res.t_hat)

    I = np.zeros_like(omega_inf)
    for i in range(1, len(t)):
        I[i] = I[i - 1] + 0.5 * (omega_inf[i] + omega_inf[i - 1]) * (t[i] - t[i - 1])

    E_shell_all = np.array(E_shell_all)
    j_thr = []
    j_pct = []
    for e_shell, total in zip(E_shell_all, total_energy):
        j_thr.append(frontier_scale(e_shell, 1e-6 * total))
        j_pct.append(percentile_frontier(e_shell, 0.05))
    j_thr = np.array(j_thr)
    j_pct = np.array(j_pct)

    j_thr_mon = monotone_envelope(j_thr)
    j_pct_mon = monotone_envelope(j_pct)

    j_mon = j_thr_mon if np.max(j_thr_mon) > np.max(j_pct_mon) else j_pct_mon
    j_vals, t_reach = crossing_times(t, j_mon)

    print("NS-BU-03 BKM proxy vs frontier")
    print(
        f"N={N} dt={dt} t_max={t_max} nu={nu} mu={mu} alpha={alpha} k_init={k_init} record_u_hat_every={record_u_hat_every}"
    )
    print(f"forcing: k_force={k_force} amp={amp}")
    print("NOTE: ||omega||_inf is grid sup; resolution-limited; de-aliasing on; diagnostic only.")
    print(f"max_omega_inf={float(np.max(omega_inf)):.6f} I_final={float(I[-1]):.6f}")
    if t_reach.size:
        print("j_new | t_cross | omega_inf | I(t_cross)")
        for j, tc in zip(j_vals, t_reach):
            idx = int(np.searchsorted(t, tc))
            print(f"{j:4d} | {tc:7.4f} | {omega_inf[idx]:9.6f} | {I[idx]:10.6f}")
    else:
        print("no crossings observed")

    fig, axes = plt.subplots(2, 1, figsize=(7, 8), sharex=True)

    axes[0].plot(t, omega_inf, label=r"$||\omega||_\infty$")
    axes[0].plot(t, I, label=r"$\int_0^t ||\omega||_\infty ds$")
    for tc in t_reach:
        axes[0].axvline(tc, color="k", linestyle=":", alpha=0.4)
    axes[0].set_ylabel("BKM proxy")
    axes[0].legend(loc="best")

    axes[1].plot(t, j_thr_mon, label="frontier thr")
    axes[1].plot(t, j_pct_mon, label="frontier pct")
    for j, tc in zip(j_vals, t_reach):
        axes[1].axvline(tc, color="k", linestyle=":", alpha=0.4)
        axes[1].text(tc, j + 0.1, f"j={j}", rotation=90, fontsize=8, va="bottom")
    axes[1].set_xlabel("t")
    axes[1].set_ylabel("j*(t)")
    axes[1].legend(loc="best")

    plt.tight_layout()
    out_path = ART / "nsbu03_bkm_frontier.png"
    plt.savefig(out_path, dpi=150)
    print(f"saved {out_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
