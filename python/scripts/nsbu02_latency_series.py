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


def compute_frontiers(u_hat_hist: np.ndarray, k_mag: np.ndarray, *, L: float) -> tuple[np.ndarray, np.ndarray, int]:
    j_max = None
    E_shell_list = []
    total_energy = []
    for u_hat in u_hat_hist:
        e_shell = shell_energies(u_hat, k_mag, j_max=j_max, L=L)
        if j_max is None:
            j_max = len(e_shell) - 1
        E_shell_list.append(e_shell)
        total_energy.append(float(np.sum(e_shell)))
    E_shell_arr = np.array(E_shell_list)
    total_energy = np.array(total_energy)

    j_thr = []
    j_pct = []
    for e_shell, total in zip(E_shell_arr, total_energy):
        thresh = 1e-6 * total
        j_thr.append(frontier_scale(e_shell, thresh))
        j_pct.append(percentile_frontier(e_shell, 0.01))
    return np.array(j_thr), np.array(j_pct), j_max


def print_latency_table(label: str, t: np.ndarray, j_series: np.ndarray, *, j_max: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    j_mon = monotone_envelope(j_series)
    j_vals, t_reach = crossing_times(t, j_mon, j_end=j_max)
    reached = ~np.isnan(t_reach)
    j_vals = j_vals[reached]
    t_reach = t_reach[reached]

    delta_t = np.diff(t_reach)
    cumulative = t_reach - t_reach[0]

    j0 = int(j_vals[0]) if j_vals.size else -1
    j_last = int(j_vals[-1]) if j_vals.size else -1
    print(f"{label}: j0={j0} j_last={j_last} j_max={j_max}")
    print("j | t_reach | delta_t | cumulative")
    for idx, j in enumerate(j_vals):
        dt_val = 0.0 if idx == 0 else float(delta_t[idx - 1])
        print(f"{j:2d} | {t_reach[idx]:7.4f} | {dt_val:7.4f} | {cumulative[idx]:9.4f}")

    if delta_t.size:
        print(f"min_delta_t={float(np.min(delta_t)):.4f} median_delta_t={float(np.median(delta_t)):.4f}")
    else:
        print("min_delta_t=nan median_delta_t=nan")

    plateau = j_last < j_max
    last_cross = t_reach[-1] if t_reach.size else 0.0
    if plateau or (t.size and last_cross < t[-1] * 0.8):
        print("WARNING: plateau or no recent crossings (finite resolution likely caps frontier)")
    print(f"monotone_steps={int(j_vals.size - 1)} last_cross={last_cross:.4f}")
    return j_mon, j_vals, t_reach


def run_case(
    *,
    N: int,
    L: float,
    dt: float,
    t_max: float,
    nu: float,
    mu: float,
    alpha: float,
    u0_hat: np.ndarray,
    forcing_hat: np.ndarray | None,
    record_u_hat_every: int,
):
    kx, ky, kz, _ = kgrid_3d(N, L=L)

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
        forcing_hat_fn=forcing_hat_fn if forcing_hat is not None else None,
        L=L,
        dealias=True,
        record_every=10,
        record_u_hat_every=record_u_hat_every,
    )
    if res.u_hat_hist is None or res.t_hat is None:
        raise RuntimeError("u_hat_hist not recorded")
    return res.t_hat, res.u_hat_hist


def main() -> int:
    N = 32
    L = 2 * np.pi
    dt = 5e-4
    t_max = 0.2
    nu_list = [0.004, 0.002, 0.001]
    mu = 0.0
    alpha = 2.0
    seed = 0
    k_init = 3
    record_u_hat_every = 20
    forcing_on = True
    k_force = 2
    amp = 1.0e-2

    ART.mkdir(parents=True, exist_ok=True)

    u0_hat = make_deterministic_u0_hat_3d(N, L=L, seed=seed, k_init=k_init, dealias=True)

    forcing_hat = None
    if forcing_on:
        forcing_hat = make_lowk_divfree_forcing_hat_3d(
            N,
            L=L,
            seed=seed,
            k_force=k_force,
            amp=amp,
            dealias=True,
        )
        forcing_label = f"forcing: k_force={k_force} amp={amp}"
    else:
        forcing_label = "forcing: none"

    kx, ky, kz, k2 = kgrid_3d(N, L=L)
    k_mag = np.sqrt(k2)

    print("NS-BU-02 latency series")
    print(f"N={N} dt={dt} t_max={t_max} nu_list={nu_list} mu={mu} alpha={alpha} record_u_hat_every={record_u_hat_every}")
    print(f"{forcing_label}")

    j_mon_pct = {}
    t_reach_pct = {}

    for nu in nu_list:
        t, u_hat_hist = run_case(
            N=N,
            L=L,
            dt=dt,
            t_max=t_max,
            nu=nu,
            mu=mu,
            alpha=alpha,
            u0_hat=u0_hat,
            forcing_hat=forcing_hat,
            record_u_hat_every=record_u_hat_every,
        )

        j_thr, j_pct, j_max = compute_frontiers(u_hat_hist, k_mag, L=L)

        print(f"\nnu={nu} (threshold frontier)")
        j_mon_thr, j_vals_thr, t_reach_thr = print_latency_table("thr", t, j_thr, j_max=j_max)
        print(f"nu={nu} (percentile frontier)")
        j_mon, j_vals_pct, t_reach = print_latency_table("pct", t, j_pct, j_max=j_max)

        j_mon_pct[nu] = j_mon
        t_reach_pct[nu] = (j_vals_pct, t_reach, j_max)

    fig, axes = plt.subplots(2, 1, figsize=(7, 8), sharex=False)

    for nu in nu_list:
        axes[0].plot(t, j_mon_pct[nu], label=f"nu={nu}")
    axes[0].set_xlabel("t")
    axes[0].set_ylabel("j*(t) percentile")
    axes[0].set_title("Monotone frontier (percentile)")
    axes[0].legend(loc="best")

    for nu in nu_list:
        j_vals, t_reach, j_max = t_reach_pct[nu]
        if t_reach.size:
            axes[1].plot(j_vals, t_reach, marker="o", label=f"nu={nu}")
    axes[1].set_xlabel("j")
    axes[1].set_ylabel("t_reach")
    axes[1].set_title("Latency curve (percentile)")
    axes[1].legend(loc="best")

    plt.tight_layout()
    out_path = ART / "nsbu02_latency.png"
    plt.savefig(out_path, dpi=150)
    print(f"saved {out_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
