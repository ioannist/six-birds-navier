import argparse
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "python" / "src"
sys.path.insert(0, str(SRC))

from nswave.ns3d_spectral import make_deterministic_u0_hat_3d, simulate_ns3d  # noqa: E402
from nswave.operators import kgrid_3d  # noqa: E402
from nswave.deterministic_npz import savez_deterministic  # noqa: E402


def energy_shell_spectrum(u_hat: np.ndarray, k2: np.ndarray) -> np.ndarray:
    k_mag = np.sqrt(k2)
    shell = np.floor(k_mag + 0.5).astype(int)
    e_mode = 0.5 * (
        np.abs(u_hat[0]) ** 2 + np.abs(u_hat[1]) ** 2 + np.abs(u_hat[2]) ** 2
    )
    kmax = int(shell.max())
    e_shell = np.bincount(shell.ravel(), weights=e_mode.ravel(), minlength=kmax + 1)
    return e_shell


def _meta_bytes(config: dict) -> np.ndarray:
    return np.frombuffer(json.dumps(config, sort_keys=True).encode("utf-8"), dtype=np.uint8)


def main() -> None:
    parser = argparse.ArgumentParser(description="NS-15 UV completion turns on (3D)")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out-dir", default=str(ROOT / "python" / "artifacts" / "paper"))
    args = parser.parse_args()

    N = 24
    L = 2 * np.pi
    dt = 5e-4
    t_max = 0.10
    nu = 1e-3
    alpha = 2.0
    mu_base = 0.0
    mu_comp = 3e-3
    dealias = True
    seed = args.seed
    k_init = 8

    record_u_hat_every = int(round(0.02 / dt))
    k_low = 4
    k0 = 7
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cfg = {
        "script": Path(__file__).name,
        "seed": seed,
        "N": N,
        "dt": dt,
        "t_max": t_max,
        "nu": nu,
        "mu_base": mu_base,
        "mu_comp": mu_comp,
        "alpha": alpha,
        "k_low": k_low,
        "k0": k0,
        "k_init": k_init,
        "out_dir": str(out_dir),
    }
    print(f"FIGURE_CONFIG: {json.dumps(cfg, sort_keys=True)}")

    u0_hat = make_deterministic_u0_hat_3d(
        N,
        L=L,
        seed=seed,
        k_init=k_init,
        dealias=dealias,
    )

    res0 = simulate_ns3d(
        N=N,
        dt=dt,
        t_max=t_max,
        nu=nu,
        mu=mu_base,
        alpha=alpha,
        u0_hat=u0_hat,
        forcing_hat_fn=None,
        L=L,
        dealias=dealias,
        record_every=10,
        record_u_hat_every=record_u_hat_every,
    )

    res1 = simulate_ns3d(
        N=N,
        dt=dt,
        t_max=t_max,
        nu=nu,
        mu=mu_comp,
        alpha=alpha,
        u0_hat=u0_hat,
        forcing_hat_fn=None,
        L=L,
        dealias=dealias,
        record_every=10,
        record_u_hat_every=record_u_hat_every,
    )

    if res0.u_hat_hist is None or res1.u_hat_hist is None:
        raise RuntimeError("u_hat history not recorded")

    kx, ky, kz, k2 = kgrid_3d(N, L=L)

    print("NS-15 UV completion turns on (3D)")
    print(
        f"N={N} dt={dt} t_max={t_max} nu={nu} alpha={alpha} "
        f"mu_base={mu_base} mu_comp={mu_comp} k_low={k_low} k0={k0}"
    )

    for j in range(len(res0.t_hat)):
        t = float(res0.t_hat[j])
        e0 = energy_shell_spectrum(res0.u_hat_hist[j], k2)
        e1 = energy_shell_spectrum(res1.u_hat_hist[j], k2)

        low0 = float(np.sum(e0[: k_low + 1]))
        low1 = float(np.sum(e1[: k_low + 1]))
        high0 = float(np.sum(e0[k0:]))
        high1 = float(np.sum(e1[k0:]))

        ratio_low = low1 / max(low0, 1e-30)
        ratio_high = high1 / max(high0, 1e-30)

        print(
            f"t={t:.3f} ratio_low={ratio_low:.3f} ratio_high={ratio_high:.3f} "
            f"high_base={high0:.3e} high_comp={high1:.3e}"
        )

    e0_final = energy_shell_spectrum(res0.u_hat_hist[-1], k2)
    e1_final = energy_shell_spectrum(res1.u_hat_hist[-1], k2)
    k = np.arange(len(e0_final))

    out_path = out_dir / "ns15_uv_completion_3d.png"

    plt.figure(figsize=(6, 4))
    plt.semilogy(k[1:], e0_final[1:], label="baseline (mu=0)")
    plt.semilogy(k[1:], e1_final[1:], label="completed (mu>0)")
    plt.xlabel("k")
    plt.ylabel("E(k)")
    plt.title(f"NS-15 UV completion (N={N}, nu={nu}, mu={mu_comp}, alpha={alpha})")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path)
    print(f"saved {out_path}")

    data_path = out_dir / "ns15_uv_completion_3d_data.npz"
    savez_deterministic(
        data_path,
        meta_json=_meta_bytes(cfg),
        t_hat=np.array(res0.t_hat),
        e0_final=e0_final,
        e1_final=e1_final,
        k=k,
    )
    print(f"saved data: {data_path}")


if __name__ == "__main__":
    main()
