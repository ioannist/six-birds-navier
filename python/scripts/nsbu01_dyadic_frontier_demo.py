import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "python" / "src"
sys.path.insert(0, str(SRC))

from nswave.dyadic import frontier_scale, percentile_frontier, shell_energies  # noqa: E402
from nswave.ns3d_spectral import make_deterministic_u0_hat_3d, simulate_ns3d  # noqa: E402
from nswave.operators import kgrid_3d  # noqa: E402


ART = ROOT / "python" / "artifacts"


def run_case(mu: float, *, N: int, dt: float, t_max: float, nu: float, alpha: float, u0_hat: np.ndarray, record_u_hat_every: int):
    res = simulate_ns3d(
        N=N,
        dt=dt,
        t_max=t_max,
        nu=nu,
        mu=mu,
        alpha=alpha,
        u0_hat=u0_hat,
        dealias=True,
        record_every=10,
        record_u_hat_every=record_u_hat_every,
    )
    if res.u_hat_hist is None or res.t_hat is None:
        raise RuntimeError("u_hat_hist not recorded")
    return res.t_hat, res.u_hat_hist


def main() -> int:
    N = 24
    L = 2 * np.pi
    dt = 5e-4
    t_max = 0.05
    nu = 1e-3
    alpha = 2.0
    mu_base = 0.0
    mu_comp = 3e-3
    record_u_hat_every = int(0.01 / dt)
    k_init = 8

    ART.mkdir(parents=True, exist_ok=True)

    u0_hat = make_deterministic_u0_hat_3d(N, L=L, seed=0, k_init=k_init, dealias=True)

    t_base, u_hat_base = run_case(
        mu_base,
        N=N,
        dt=dt,
        t_max=t_max,
        nu=nu,
        alpha=alpha,
        u0_hat=u0_hat,
        record_u_hat_every=record_u_hat_every,
    )
    t_comp, u_hat_comp = run_case(
        mu_comp,
        N=N,
        dt=dt,
        t_max=t_max,
        nu=nu,
        alpha=alpha,
        u0_hat=u0_hat,
        record_u_hat_every=record_u_hat_every,
    )

    kx, ky, kz, k2 = kgrid_3d(N, L=L)
    k_mag = np.sqrt(k2)

    E_base = []
    E_comp = []
    front_base = []
    front_comp = []
    pct_base = []
    pct_comp = []

    j_max = None
    for u_hat in u_hat_base:
        e_shell = shell_energies(u_hat, k_mag, j_max=j_max, L=L)
        if j_max is None:
            j_max = len(e_shell) - 1
        E_base.append(e_shell)
        thresh = 1e-6 * float(np.sum(e_shell))
        front_base.append(frontier_scale(e_shell, thresh))
        pct_base.append(percentile_frontier(e_shell, 0.05))

    for u_hat in u_hat_comp:
        e_shell = shell_energies(u_hat, k_mag, j_max=j_max, L=L)
        E_comp.append(e_shell)
        thresh = 1e-6 * float(np.sum(e_shell))
        front_comp.append(frontier_scale(e_shell, thresh))
        pct_comp.append(percentile_frontier(e_shell, 0.05))

    E_base = np.array(E_base)
    E_comp = np.array(E_comp)
    front_base = np.array(front_base)
    front_comp = np.array(front_comp)
    pct_base = np.array(pct_base)
    pct_comp = np.array(pct_comp)

    print("NS-BU-01 dyadic frontier demo")
    print(f"N={N} dt={dt} t_max={t_max} nu={nu} mu={mu_comp} alpha={alpha} record_u_hat_every={record_u_hat_every}")
    print(f"j_max={j_max}")
    print("threshold=1e-6 * total_energy")
    print(f"frontier_final_base={front_base[-1]} frontier_final_comp={front_comp[-1]}")
    print(f"percentile_final_base={pct_base[-1]} percentile_final_comp={pct_comp[-1]}")

    eps = 1e-16
    fig, axes = plt.subplots(3, 1, figsize=(7, 8), sharex=True)

    im0 = axes[0].imshow(
        np.log10(E_base + eps).T,
        aspect="auto",
        origin="lower",
        extent=[t_base[0], t_base[-1], 0, j_max],
    )
    axes[0].set_title("Baseline (mu=0): log10 E_j(t)")
    axes[0].set_ylabel("j")
    fig.colorbar(im0, ax=axes[0], fraction=0.025, pad=0.02)

    im1 = axes[1].imshow(
        np.log10(E_comp + eps).T,
        aspect="auto",
        origin="lower",
        extent=[t_comp[0], t_comp[-1], 0, j_max],
    )
    axes[1].set_title("Completed (mu>0): log10 E_j(t)")
    axes[1].set_ylabel("j")
    fig.colorbar(im1, ax=axes[1], fraction=0.025, pad=0.02)

    axes[2].plot(t_base, front_base, label="frontier (mu=0)")
    axes[2].plot(t_comp, front_comp, label="frontier (mu>0)")
    axes[2].plot(t_base, pct_base, "--", label="pct frontier (mu=0)")
    axes[2].plot(t_comp, pct_comp, "--", label="pct frontier (mu>0)")
    axes[2].set_xlabel("t")
    axes[2].set_ylabel("j*")
    axes[2].legend(loc="best")

    plt.tight_layout()
    out_path = ART / "nsbu01_frontier.png"
    plt.savefig(out_path, dpi=150)
    print(f"saved {out_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
