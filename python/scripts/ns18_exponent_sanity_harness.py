import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "python" / "src"
sys.path.insert(0, str(SRC))

from nswave.operators import kgrid_3d, project_div_free_3d  # noqa: E402
from nswave.deterministic_npz import savez_deterministic  # noqa: E402


ART_DIR = ROOT / "python" / "artifacts"


def make_shell_field_u_hat(
    N: int,
    K: int,
    *,
    L: float,
    seed: int,
    shell_halfwidth: float = 0.5,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    u = rng.standard_normal((3, N, N, N))
    u_hat = np.fft.fftn(u, axes=(1, 2, 3))
    u_hat[:, 0, 0, 0] = 0.0

    kx, ky, kz, k2 = kgrid_3d(N, L=L)
    kk = np.fft.fftfreq(N) * N
    kx_i, ky_i, kz_i = np.meshgrid(kk, kk, kk, indexing="ij")
    kmag = np.sqrt(kx_i * kx_i + ky_i * ky_i + kz_i * kz_i)
    shell_mask = np.abs(kmag - K) <= shell_halfwidth

    u_hat = u_hat * shell_mask
    u_hat = project_div_free_3d(u_hat, kx, ky, kz)
    return u_hat, kx, ky, kz, k2


def norm_lambda(u_hat: np.ndarray, k2: np.ndarray, r: float, L: float) -> float:
    mult = np.zeros_like(k2)
    mask = k2 > 0
    mult[mask] = k2[mask] ** (r / 2.0)
    ur0 = np.fft.ifftn(mult * u_hat[0]).real
    ur1 = np.fft.ifftn(mult * u_hat[1]).real
    ur2 = np.fft.ifftn(mult * u_hat[2]).real
    val = np.mean(ur0 * ur0 + ur1 * ur1 + ur2 * ur2) * (L**3)
    return float(math.sqrt(max(val, 0.0)))


def slope_loglog(x: np.ndarray, y: np.ndarray, tiny: float = 1e-30) -> float:
    xlog = np.log(np.maximum(x, tiny))
    ylog = np.log(np.maximum(y, tiny))
    coeffs = np.polyfit(xlog, ylog, 1)
    return float(coeffs[0])


def _meta_bytes(config: dict) -> np.ndarray:
    return np.frombuffer(json.dumps(config, sort_keys=True).encode("utf-8"), dtype=np.uint8)


def main() -> int:
    parser = argparse.ArgumentParser(description="NS-18 exponent sanity harness")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out-dir", default=str(ROOT / "python" / "artifacts" / "paper"))
    args = parser.parse_args()

    N = 32
    L = 2 * np.pi
    s = 3.0
    alpha = 1.25
    K_list = [4, 6, 8, 10]
    reps = 3
    eps = 0.1
    theta_good = 1.0 / alpha
    theta_bad = 0.4
    seed = args.seed
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cfg = {
        "script": Path(__file__).name,
        "seed": seed,
        "N": N,
        "dt": 0.0,
        "t_max": 0.0,
        "nu": 0.0,
        "mu": 0.0,
        "alpha": alpha,
        "s": s,
        "K_list": K_list,
        "reps": reps,
        "eps": eps,
        "theta_good": theta_good,
        "theta_bad": theta_bad,
        "out_dir": str(out_dir),
    }
    print(f"FIGURE_CONFIG: {json.dumps(cfg, sort_keys=True)}")

    def compute_for_theta(theta: float) -> tuple[np.ndarray, np.ndarray]:
        interp_vals = []
        creq_vals = []
        m = 3.0 - theta
        q = 2.0 / (2.0 - theta)

        for K in K_list:
            r_interp = []
            r_creq = []
            for rep in range(reps):
                u_hat, _kx, _ky, _kz, k2 = make_shell_field_u_hat(
                    N, K, L=L, seed=seed + rep
                )
                A = norm_lambda(u_hat, k2, s, L)
                A1 = norm_lambda(u_hat, k2, s + 1.0, L)
                B = norm_lambda(u_hat, k2, s + alpha, L)

                denom = (A ** (1.0 - theta)) * (B ** theta)
                r_interp.append(A1 / max(denom, 1e-30))

                term = (A**m) * (B**theta)
                creq = max(0.0, term - eps * (B**2)) / (A ** (m * q) + 1e-30)
                r_creq.append(creq)

            interp_vals.append(float(np.median(r_interp)))
            creq_vals.append(float(np.median(r_creq)))

        return np.array(interp_vals), np.array(creq_vals)

    interp_good, creq_good = compute_for_theta(theta_good)
    interp_bad, creq_bad = compute_for_theta(theta_bad)

    K_arr = np.array(K_list, dtype=float)
    slope_interp_good = slope_loglog(K_arr, interp_good)
    slope_interp_bad = slope_loglog(K_arr, interp_bad)
    slope_creq_good = slope_loglog(K_arr, creq_good + 1e-30)
    slope_creq_bad = slope_loglog(K_arr, creq_bad + 1e-30)

    interp_status_good = "PASS" if abs(slope_interp_good) <= 0.15 else "FAIL"
    interp_status_bad = "FAIL" if abs(slope_interp_bad) >= 0.30 else "PASS"
    ordering_ok = slope_creq_good <= slope_creq_bad + 0.1

    print("NS-18 exponent sanity harness")
    print(
        "N={N} s={s} alpha={alpha} theta_good={tg:.6f} theta_bad={tb:.6f} eps={eps}".format(
            N=N, s=s, alpha=alpha, tg=theta_good, tb=theta_bad, eps=eps
        )
    )
    print("Interpolation slopes:")
    print(
        "  theta={:.6f} slope_interp={:.3f} {}".format(
            theta_good, slope_interp_good, interp_status_good
        )
    )
    print(
        "  theta={:.6f} slope_interp={:.3f} {}".format(
            theta_bad, slope_interp_bad, interp_status_bad
        )
    )
    print("Young-absorption proxy slopes:")
    print("  theta={:.6f} slope_creq={:.3f}".format(theta_good, slope_creq_good))
    print("  theta={:.6f} slope_creq={:.3f}".format(theta_bad, slope_creq_bad))
    print(f"  ordering_ok={ordering_ok}")
    print("Table (K, interp_good, interp_bad)")
    for K, ig, ib in zip(K_list, interp_good, interp_bad):
        print(f"  {K:>2d} {ig:.6e} {ib:.6e}")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].loglog(K_arr, interp_good, "o-", label=f"theta={theta_good:.3f}")
    axes[0].loglog(K_arr, interp_bad, "o-", label=f"theta={theta_bad:.3f}")
    axes[0].set_xlabel("K")
    axes[0].set_ylabel("R_interp")
    axes[0].set_title("Interpolation ratio")
    axes[0].legend()

    axes[1].loglog(K_arr, creq_good + 1e-30, "o-", label=f"theta={theta_good:.3f}")
    axes[1].loglog(K_arr, creq_bad + 1e-30, "o-", label=f"theta={theta_bad:.3f}")
    axes[1].set_xlabel("K")
    axes[1].set_ylabel("C_req")
    axes[1].set_title("Young proxy")
    axes[1].legend()

    fig.suptitle(f"NS-18 exponent sanity (N={N}, alpha={alpha})")
    fig.tight_layout()

    out_path = out_dir / "ns18_exponent_sanity.png"
    fig.savefig(out_path, dpi=150)
    print(f"saved {out_path}")

    data_path = out_dir / "ns18_exponent_sanity_data.npz"
    savez_deterministic(
        data_path,
        meta_json=_meta_bytes(cfg),
        K_arr=K_arr,
        interp_good=interp_good,
        interp_bad=interp_bad,
        creq_good=creq_good,
        creq_bad=creq_bad,
        slope_interp_good=np.array(slope_interp_good),
        slope_interp_bad=np.array(slope_interp_bad),
        slope_creq_good=np.array(slope_creq_good),
        slope_creq_bad=np.array(slope_creq_bad),
        ordering_ok=np.array(ordering_ok),
    )
    print(f"saved data: {data_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
