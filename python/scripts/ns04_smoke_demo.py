import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "python" / "src"
sys.path.insert(0, str(SRC))

from nswave.ns2d_spectral import simulate_ns2d  # noqa: E402


def main() -> None:
    res = simulate_ns2d(
        N=32,
        dt=5e-4,
        t_max=0.05,
        nu=1e-3,
        mu=1e-4,
        alpha=2.0,
        forcing_hat_fn=None,
        record_every=10,
    )

    print("idx   t       enstrophy    diss_nu_E   diss_mu_E")
    for i in range(len(res.t)):
        print(
            f"{i:3d}  {res.t[i]:.4f}  {res.enstrophy[i]:.6f}  "
            f"{res.diss_nu[i]:.6f}  {res.diss_mu[i]:.6f}"
        )


if __name__ == "__main__":
    main()
