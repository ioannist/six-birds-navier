import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "python" / "src"
sys.path.insert(0, str(SRC))

from nswave.ns3d_spectral import make_deterministic_u0_hat_3d, simulate_ns3d  # noqa: E402


def main() -> None:
    N = 16
    dt = 5e-4
    t_max = 0.02
    nu = 1e-3
    mu = 1e-4
    alpha = 2.0

    u0_hat = make_deterministic_u0_hat_3d(N, seed=0, k_init=6)
    res = simulate_ns3d(
        N=N,
        dt=dt,
        t_max=t_max,
        nu=nu,
        mu=mu,
        alpha=alpha,
        u0_hat=u0_hat,
        record_every=10,
    )

    print("idx  t       energy      diss_nu    diss_mu    div_rms")
    for i in range(len(res.t)):
        print(
            f"{i:3d}  {res.t[i]:.4f}  {res.energy[i]:.6e}  "
            f"{res.diss_nu[i]:.6e}  {res.diss_mu[i]:.6e}  {res.div_rms[i]:.3e}"
        )

    print(
        f"E0={res.energy[0]:.6e} Efinal={res.energy[-1]:.6e} "
        f"max_div_rms={res.div_rms.max():.3e} success={res.success} message={res.message}"
    )


if __name__ == "__main__":
    main()
