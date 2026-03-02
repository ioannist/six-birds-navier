import argparse
import json
import sys
from pathlib import Path
import subprocess

import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "python" / "src"
sys.path.insert(0, str(SRC))

from nswave.deterministic_npz import savez_deterministic  # noqa: E402
from nswave.positive_real_fit import debye_model, fit_debye_positive_real  # noqa: E402

ART = ROOT / "python" / "artifacts"


def load_data(data_path: Path, *, seed: int, out_dir: Path) -> dict:
    if not data_path.exists():
        print("nsbu04 data not found; running nsbu04 to generate it")
        proc = subprocess.run(
            [
                sys.executable,
                "python/scripts/nsbu04_sgs_frequency_response.py",
                "--seed",
                str(seed),
                "--out-dir",
                str(out_dir),
            ],
            cwd=ROOT,
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            print(proc.stdout)
            print(proc.stderr)
            raise RuntimeError("nsbu04 failed; cannot proceed")
    return dict(np.load(data_path, allow_pickle=True))


def main() -> int:
    parser = argparse.ArgumentParser(description="NS-BU-05 SGS realization fit")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out-dir", default=str(ROOT / "python" / "artifacts" / "paper"))
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    data_path = out_dir / "nsbu04_sgs_H_data.npz"
    data = load_data(data_path, seed=args.seed, out_dir=out_dir)

    omega = data["omega"]
    H = data["H"]
    denom = data["denom"]
    shells = data["shells"]

    cfg = {
        "script": Path(__file__).name,
        "seed": args.seed,
        "N": 0,
        "dt": 0.0,
        "t_max": 0.0,
        "nu": 0.0,
        "mu": 0.0,
        "alpha": 0.0,
        "source_npz": str(data_path),
        "out_dir": str(out_dir),
    }
    print(f"FIGURE_CONFIG: {json.dumps(cfg, sort_keys=True)}")

    print("NS-BU-05 SGS realization fit")
    print(f"loaded {data_path}")

    denom_tol = 1e-12
    passivity_rows = []
    selected_shells = []
    for i, s in enumerate(shells):
        valid = (denom[i] > denom_tol) & (omega > 0)
        if np.any(valid):
            re_vals = np.real(H[i][valid])
            re_med = float(np.median(re_vals))
            frac_nonneg = float(np.mean(re_vals >= -1e-3))
            passivity_rows.append((int(s), re_med, frac_nonneg))
            if frac_nonneg >= 0.6 and re_med > 0:
                selected_shells.append(i)
        else:
            passivity_rows.append((int(s), float("nan"), float("nan")))

    print("shell | re_median | frac_nonneg")
    for s, re_med, frac_nonneg in passivity_rows:
        print(f"{s:5d} | {re_med:9.4f} | {frac_nonneg:11.4f}")

    if not selected_shells:
        print("NO SHELLS PASSED SELECTION")
        return 0

    fit_results = {}
    for idx in selected_shells:
        s = int(shells[idx])
        valid = (denom[idx] > denom_tol) & (omega > 0)
        omega_pos = omega[valid]
        H_pos = H[idx][valid]
        weights = denom[idx][valid]
        for n_modes in [1, 2, 3, 4]:
            res = fit_debye_positive_real(
                omega_pos,
                H_pos,
                n_modes=n_modes,
                weights=weights,
                max_nfev=500,
                seed=0,
            )
            fit_results[(s, n_modes)] = res

    print("shell | M | rmse | rmse_w | min_Re_fit | frac_nonneg_fit")
    for (s, n_modes), res in fit_results.items():
        print(
            f"{s:5d} | {n_modes:1d} | {res.rmse:7.4f} | {res.rmse_weighted:7.4f} | "
            f"{res.min_re_fit:11.4f} | {res.frac_re_nonneg_fit:15.4f}"
        )

    # choose representative shell: lowest selected shell
    rep_idx = selected_shells[0]
    rep_shell = int(shells[rep_idx])
    best = None
    for n_modes in [1, 2, 3, 4]:
        res = fit_results[(rep_shell, n_modes)]
        if best is None:
            best = res
            best_M = n_modes
            continue
        if res.rmse_weighted < best.rmse_weighted - 1e-12:
            best = res
            best_M = n_modes
        elif abs(res.rmse_weighted - best.rmse_weighted) <= 1e-12 and n_modes < best_M:
            best = res
            best_M = n_modes
            print("best_M chosen with tie-break: smallest M")

    valid = (denom[rep_idx] > denom_tol) & (omega > 0)
    omega_pos = omega[valid]
    H_pos = H[rep_idx][valid]

    H_fit = debye_model(omega_pos, best.a0, best.a, best.b)
    omega_dense = np.logspace(np.log10(np.min(omega_pos)), np.log10(np.max(omega_pos)), 400)
    H_dense = debye_model(omega_dense, best.a0, best.a, best.b)
    min_re_dense = float(np.min(np.real(H_dense)))
    if min_re_dense < -1e-6:
        raise RuntimeError("fit is not positive-real on dense grid")

    print(f"selected_shell={rep_shell} best_M={best_M}")
    print(f"best_params: a0={best.a0:.6f} a={best.a} b={best.b}")

    fig, axes = plt.subplots(3, 1, figsize=(7, 8), sharex=False)
    axes[0].plot(omega_pos, np.real(H_pos), label="Re H data")
    axes[0].plot(omega_pos, np.real(H_fit), label="Re H fit")
    axes[0].axhline(0.0, color="k", linestyle="--", linewidth=0.8)
    axes[0].set_ylabel("Re H")
    axes[0].legend(loc="best")

    axes[1].plot(omega_pos, np.abs(H_pos - H_fit))
    axes[1].set_ylabel("|H-H_fit|")

    rmse_vals = [fit_results[(rep_shell, m)].rmse_weighted for m in [1, 2, 3, 4]]
    axes[2].plot([1, 2, 3, 4], rmse_vals, marker="o")
    axes[2].set_xlabel("# modes")
    axes[2].set_ylabel("weighted RMSE")

    plt.tight_layout()
    out_path = out_dir / "nsbu05_sgs_fit.png"
    plt.savefig(out_path, dpi=150)
    print(f"saved {out_path}")

    data_out = out_dir / "nsbu05_sgs_fit_data.npz"
    meta_json = np.frombuffer(json.dumps(cfg, sort_keys=True).encode("utf-8"), dtype=np.uint8)
    savez_deterministic(
        data_out,
        meta_json=meta_json,
        omega_pos=omega_pos,
        H_pos=H_pos,
        H_fit=H_fit,
        omega_dense=omega_dense,
        H_dense=H_dense,
        rmse_vals=np.array(rmse_vals),
        best_M=np.array(best_M),
        best_a0=np.array(best.a0),
        best_a=best.a,
        best_b=best.b,
    )
    print(f"saved data: {data_out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
