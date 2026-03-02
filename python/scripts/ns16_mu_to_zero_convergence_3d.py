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


def low_k_error(u_hat_base, u_hat_mu, k2, L, k_low):
    k_mag = np.sqrt(k2)
    mask = k_mag <= k_low

    u0 = []
    u1 = []
    for i in range(3):
        u0.append(np.fft.ifftn(u_hat_base[i] * mask).real)
        u1.append(np.fft.ifftn(u_hat_mu[i] * mask).real)

    norm_base = np.sqrt(np.mean(u0[0] ** 2 + u0[1] ** 2 + u0[2] ** 2) * L**3)
    diff = (u1[0] - u0[0]) ** 2 + (u1[1] - u0[1]) ** 2 + (u1[2] - u0[2]) ** 2
    norm_err = np.sqrt(np.mean(diff) * L**3)
    return float(norm_err / max(norm_base, 1e-12))


def _meta_bytes(config: dict) -> np.ndarray:
    return np.frombuffer(json.dumps(config, sort_keys=True).encode("utf-8"), dtype=np.uint8)


def main() -> None:
    parser = argparse.ArgumentParser(description="NS-16 mu->0 convergence (3D)")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out-dir", default=str(ROOT / "python" / "artifacts" / "paper"))
    args = parser.parse_args()

    N = 24
    L = 2 * np.pi
    dt = 5e-4
    t_max = 0.10
    nu = 1e-3
    alpha = 2.0
    mu_list = [0.0, 1e-2, 3e-3, 1e-3]
    dealias = True
    seed = args.seed
    k_init = 8
    k_low = 4
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cfg = {
        "script": Path(__file__).name,
        "seed": seed,
        "N": N,
        "dt": dt,
        "t_max": t_max,
        "nu": nu,
        "mu": 0.0,
        "alpha": alpha,
        "mu_list": mu_list,
        "k_init": k_init,
        "k_low": k_low,
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
    kx, ky, kz, k2 = kgrid_3d(N, L=L)

    res0 = simulate_ns3d(
        N=N,
        dt=dt,
        t_max=t_max,
        nu=nu,
        mu=0.0,
        alpha=alpha,
        u0_hat=u0_hat,
        forcing_hat_fn=None,
        L=L,
        dealias=dealias,
        record_every=10,
    )
    u_hat_base = res0.u_hat_final

    errs = {}
    for mu in mu_list[1:]:
        res = simulate_ns3d(
            N=N,
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
        )
        errs[mu] = low_k_error(u_hat_base, res.u_hat_final, k2, L, k_low)

    print("NS-16 mu->0 convergence (3D low-k)")
    print(
        f"N={N} dt={dt} t_max={t_max} nu={nu} alpha={alpha} k_low={k_low} k_init={k_init}"
    )
    print("mu | err_lowk")
    print("0.0 | 0.0")
    for mu in mu_list[1:]:
        print(f"{mu:.1e} | {errs[mu]:.6e}")

    mu_sorted = sorted(errs.keys(), reverse=True)
    monotone = all(errs[mu_sorted[i + 1]] <= errs[mu_sorted[i]] + 1e-12 for i in range(len(mu_sorted) - 1))
    print(f"monotone_decreasing = {monotone}")

    out_path = out_dir / "ns16_mu_convergence_3d.png"

    mu_vals = np.array(mu_sorted)[::-1]
    err_vals = np.array([errs[mu] for mu in mu_vals])

    plt.figure(figsize=(5, 4))
    plt.loglog(mu_vals, err_vals, marker="o")
    plt.xlabel("mu")
    plt.ylabel("low-k relative error")
    plt.title(f"NS-16 mu->0 convergence (N={N}, nu={nu}, alpha={alpha}, k_low={k_low})")
    plt.tight_layout()
    plt.savefig(out_path)
    print(f"saved {out_path}")

    data_path = out_dir / "ns16_mu_convergence_3d_data.npz"
    savez_deterministic(
        data_path,
        meta_json=_meta_bytes(cfg),
        mu_vals=mu_vals,
        err_vals=err_vals,
        mu_sorted=np.array(mu_sorted),
        monotone_decreasing=np.array(monotone),
    )
    print(f"saved data: {data_path}")


if __name__ == "__main__":
    main()
