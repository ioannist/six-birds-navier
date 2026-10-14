import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "python" / "src"
sys.path.insert(0, str(SRC))

from nswave.galerkin_core import (  # noqa: E402
    B_of_u,
    GalerkinModel,
    dissipation,
    energy,
    make_energy_preserving_tensor,
    make_psd_operator,
    step_rk4,
)


def main() -> None:
    n = 8
    rng = np.random.default_rng(0)
    T = rng.standard_normal((n, n, n))
    C = make_energy_preserving_tensor(T)
    A = make_psd_operator(n, seed=1)
    Aalpha = make_psd_operator(n, seed=2)

    model = GalerkinModel(n=n, A=A, Aalpha=Aalpha, nu=0.5, mu=0.2, C=C)
    u = rng.standard_normal(n)

    dt = 1e-3
    nsteps = 100

    print("step   E        dissipation   |u·B|")
    for step in range(nsteps + 1):
        if step % 10 == 0:
            Bu = B_of_u(model, u)
            val = abs(float(u @ Bu))
            print(f"{step:4d}  {energy(u):.6f}  {dissipation(model, u):.6f}     {val:.3e}")
        if step < nsteps:
            u = step_rk4(model, u, dt)


if __name__ == "__main__":
    main()
