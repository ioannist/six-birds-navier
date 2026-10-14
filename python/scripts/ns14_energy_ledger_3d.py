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
from nswave.deterministic_npz import savez_deterministic  # noqa: E402


def cumulative_trapz(t: np.ndarray, y: np.ndarray) -> np.ndarray:
    out = np.zeros_like(y)
    if len(t) < 2:
        return out
    dt = t[1:] - t[:-1]
    out[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * dt)
    return out


def _meta_bytes(config: dict) -> np.ndarray:
    return np.frombuffer(json.dumps(config, sort_keys=True).encode("utf-8"), dtype=np.uint8)


def main() -> None:
    parser = argparse.ArgumentParser(description="NS-14 energy ledger (3D)")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out-dir", default=str(ROOT / "python" / "artifacts" / "paper"))
    args = parser.parse_args()

    N = 16
    L = 2 * np.pi
    dt = 5e-4
    t_max = 0.05
    nu = 1e-3
    mu = 1e-4
    alpha = 2.0
    seed = args.seed
    k_init = 6
    dealias = True
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cfg = {
        "script": Path(__file__).name,
        "seed": seed,
        "N": N,
        "dt": dt,
        "t_max": t_max,
        "nu": nu,
        "mu": mu,
        "alpha": alpha,
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
        record_every=1,
    )

    t = res.t
    E = res.energy
    D = res.diss_nu + res.diss_mu
    I = cumulative_trapz(t, D)

    E0 = E[0]
    denom = max(E0, 1e-12)
    R = (E + I - E0) / denom

    resid_final = float(R[-1])
    max_abs_resid = float(np.max(np.abs(R)))
    t0 = 0.005
    mask = t >= t0
    max_abs_resid_after_t0 = float(np.max(np.abs(R[mask]))) if np.any(mask) else max_abs_resid

    print("NS-14 energy ledger (3D, unforced)")
    print(f"N={N} dt={dt} t_max={t_max} nu={nu} mu={mu} alpha={alpha} k_init={k_init}")
    print(f"E0={E0:.6e} Efinal={E[-1]:.6e}")
    print(f"I_final={I[-1]:.6e}")
    print(f"resid_final={resid_final:.6e}")
    print(f"max_abs_resid={max_abs_resid:.6e}")
    print(f"max_abs_resid_after_t0={max_abs_resid_after_t0:.6e}")

    out_path = out_dir / "ns14_energy_ledger_3d.png"

    plt.figure(figsize=(6, 4))
    plt.plot(t, E / denom, label="E(t)/E0")
    plt.plot(t, R, label="ledger residual")
    plt.xlabel("t")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path)
    print(f"saved {out_path}")

    data_path = out_dir / "ns14_energy_ledger_3d_data.npz"
    savez_deterministic(
        data_path,
        meta_json=_meta_bytes(cfg),
        t=t,
        E=E,
        D=D,
        I=I,
        R=R,
        E0=np.array(E0),
        resid_final=np.array(resid_final),
        max_abs_resid=np.array(max_abs_resid),
        max_abs_resid_after_t0=np.array(max_abs_resid_after_t0),
    )
    print(f"saved data: {data_path}")


if __name__ == "__main__":
    main()
