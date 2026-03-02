import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "python" / "src"
sys.path.insert(0, str(SRC))

from nswave.ns2d_spectral import simulate_ns2d  # noqa: E402


def cumulative_trapz(t: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Cumulative trapezoid integral with y, t arrays."""
    out = np.zeros_like(y)
    if len(t) < 2:
        return out
    dt = t[1:] - t[:-1]
    out[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * dt)
    return out


def main() -> None:
    N = 32
    dt = 5e-4
    t_max = 0.2
    nu = 1e-3
    mu = 1e-4
    alpha = 2.0

    res = simulate_ns2d(
        N=N,
        dt=dt,
        t_max=t_max,
        nu=nu,
        mu=mu,
        alpha=alpha,
        forcing_hat_fn=None,
        dealias=True,
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

    t0 = 0.02
    mask = t >= t0
    max_abs_resid_after_t0 = float(np.max(np.abs(R[mask]))) if np.any(mask) else max_abs_resid

    print("NS-05 energy ledger (unforced)")
    print(f"N={N} dt={dt} t_max={t_max} nu={nu} mu={mu} alpha={alpha}")
    print(f"E0={E0:.6e} Efinal={E[-1]:.6e}")
    print(f"I_final={I[-1]:.6e}")
    print(f"resid_final={resid_final:.6e}")
    print(f"max_abs_resid={max_abs_resid:.6e}")
    print(f"max_abs_resid_after_t0={max_abs_resid_after_t0:.6e}")

    artifacts = ROOT / "python" / "artifacts"
    artifacts.mkdir(parents=True, exist_ok=True)
    out_path = artifacts / "ns05_energy_ledger.png"

    plt.figure(figsize=(6, 4))
    plt.plot(t, E / denom, label="E(t)/E0")
    plt.plot(t, R, label="ledger residual")
    plt.xlabel("t")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path)
    print(f"saved {out_path}")


if __name__ == "__main__":
    main()
