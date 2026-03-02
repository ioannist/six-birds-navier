import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "python" / "src"
sys.path.insert(0, str(SRC))

from nswave.deterministic_npz import savez_deterministic
from nswave.route_capacity_hypothesis import fit_smallest_p, assess_mismatch_summability


ART = ROOT / "python" / "artifacts"


def main() -> int:
    parser = argparse.ArgumentParser(description="NS-BU-23 route-capacity hypothesis check")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out-dir", default=str(ROOT / "python" / "artifacts" / "paper"))
    parser.add_argument("--npz", default=None)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    data_path = Path(args.npz) if args.npz is not None else out_dir / "nsbu21_route_capacity_numbers_data.npz"
    cfg = {
        "script": Path(__file__).name,
        "seed": args.seed,
        "N": 0,
        "dt": 0.0,
        "t_max": 0.0,
        "nu": 0.0,
        "mu": 0.0,
        "alpha": 0.0,
        "npz": str(data_path),
        "out_dir": str(out_dir),
    }
    print(f"FIGURE_CONFIG: {json.dumps(cfg, sort_keys=True)}")

    if not data_path.exists():
        print("data not found; running nsbu21_route_capacity_numbers.py to generate...")
        script = ROOT / "python" / "scripts" / "nsbu21_route_capacity_numbers.py"
        result = subprocess.run(
            [
                sys.executable,
                str(script),
                "--seed",
                str(args.seed),
                "--out-dir",
                str(out_dir),
            ],
            check=False,
        )
        if result.returncode != 0:
            return result.returncode

    data = np.load(data_path)
    c_tilde = data["C_tilde"]
    e = data["e_route"] if "e_route" in data else data["e"]
    e_tilde = data["e_route_tilde"] if "e_route_tilde" in data else data["e_tilde"]
    slack_poswork = data["slack_poswork"] if "slack_poswork" in data else None
    partial_sum_e = data["partial_sum_e_route"] if "partial_sum_e_route" in data else data["partial_sum_e"]
    u_rms = float(data["U_rms"])
    k_cut = data["k_cut"]
    if e_tilde is None:
        e_tilde = e / np.maximum(u_rms * k_cut, 1e-12)

    eps = 1e-12
    c_tol = 1e-14
    start_idx = int(np.argmax(c_tilde > c_tol))
    if start_idx >= len(c_tilde) - 1:
        print("insufficient c_tilde entries above tolerance; cannot compute r_steps")
        return 1

    r_steps = c_tilde[start_idx + 1 :] / np.maximum(c_tilde[start_idx:-1], eps)
    p_best, max_violation = fit_smallest_p(
        r_steps,
        margin_rel=0.05,
        p_max=10,
        step_offset=start_idx,
    )

    if p_best is not None:
        steps = np.arange(len(r_steps)) + start_idx
        r_model = ((steps + 2) / (steps + 1)) ** p_best
        worst_step = int(np.argmax(r_steps / np.maximum(r_model, eps)))
    else:
        worst_step = int(np.argmax(r_steps))
    e_tilde = e / np.maximum(u_rms * k_cut, eps)
    psum_e_tilde = np.cumsum(e_tilde)

    psum_e_final = float(partial_sum_e[-1])
    psum_e_tilde_final = float(psum_e_tilde[-1])

    summability = assess_mismatch_summability(
        e,
        e_tilde,
        data["Lambda_direct"] if "Lambda_direct" in data else np.zeros_like(e),
        tol_abs=1e-14,
        tol_rel=1e-8,
        min_steps=8,
    )
    summability_status = summability["status"]

    status_pass = summability_status.startswith("PASS")
    status_fail = summability_status == "FAIL"
    verdict = "INCONCLUSIVE"
    if p_best is not None and status_pass:
        verdict = "PLAUSIBLE"
    elif p_best is None and status_fail:
        verdict = "NOT PLAUSIBLE"

    print("NS-BU-23 route-capacity hypothesis check")
    print(f"start_idx={start_idx} c_tilde_start={c_tilde[start_idx]:.6e}")
    print(f"r_steps_first={r_steps[:3]}")
    margin_rel = 0.05
    print(f"p_best={p_best} max_violation={max_violation:.6e} margin={margin_rel:.3f} worst_step={worst_step}")
    print(f"r_steps_summary: max={float(np.max(r_steps)):.6e} median={float(np.median(r_steps)):.6e}")
    print(f"psum_e_final={psum_e_final:.6e} psum_e_tilde_final={psum_e_tilde_final:.6e}")
    print(f"MISMATCH_SUMMABILITY: {summability_status}")
    print(
        f"mismatch_tol_abs={summability['tol_abs']:.3e} mismatch_tol_rel={summability['tol_rel']:.3e} "
        f"mismatch_scale={summability['scale']:.3e}"
    )
    q_raw = "SKIP" if summability["q_est_raw"] is None else f"{summability['q_est_raw']:.3f}"
    q_tilde = "SKIP" if summability["q_est_tilde"] is None else f"{summability['q_est_tilde']:.3f}"
    tr_raw = "SKIP" if summability["tail_ratio_median_raw"] is None else f"{summability['tail_ratio_median_raw']:.3f}"
    tr_tilde = (
        "SKIP"
        if summability["tail_ratio_median_tilde"] is None
        else f"{summability['tail_ratio_median_tilde']:.3f}"
    )
    print(f"K={summability['K']} K_tail={summability['K_tail']} q_est_raw={q_raw} q_est_tilde={q_tilde}")
    print(f"tail_ratio_median_raw={tr_raw} tail_ratio_median_tilde={tr_tilde}")
    print("TAME_STEP_FACTOR: " + ("PASS" if p_best is not None else "FAIL"))
    if slack_poswork is not None:
        print(
            f"POSWORK_SLACK: min={float(np.min(slack_poswork)):.3e} "
            f"max={float(np.max(slack_poswork)):.3e}"
        )
    print(f"OVERALL: {verdict}")

    fig, axes = plt.subplots(3, 1, figsize=(7, 9))
    steps = np.arange(len(r_steps)) + start_idx
    axes[0].plot(steps, r_steps, marker="o", label="r_tilde")
    if p_best is not None:
        r_model = ((steps + 2) / (steps + 1)) ** p_best
        axes[0].plot(steps, r_model, linestyle="--", label=f"p={p_best}")
        if p_best - 1 >= 0:
            axes[0].plot(
                steps,
                ((steps + 2) / (steps + 1)) ** (p_best - 1),
                linestyle="--",
                label=f"p={p_best-1}",
            )
        axes[0].plot(
            steps,
            ((steps + 2) / (steps + 1)) ** (p_best + 1),
            linestyle="--",
            label=f"p={p_best+1}",
        )
    axes[0].set_xlabel("step index")
    axes[0].set_ylabel("r_tilde")
    axes[0].legend(loc="best")

    axes[1].plot(np.arange(len(e)), e, marker="o", label="e_route")
    if slack_poswork is not None:
        axes[1].plot(np.arange(len(slack_poswork)), slack_poswork, marker="s", label="slack_poswork")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("depth index")
    axes[1].set_ylabel("mismatch/slack (log)")
    axes[1].legend(loc="best")

    axes[2].plot(np.arange(len(partial_sum_e)), partial_sum_e, marker="o", label="sum e_route")
    axes[2].plot(np.arange(len(psum_e_tilde)), psum_e_tilde, marker="s", label="sum e_route_tilde")
    axes[2].set_xlabel("depth index")
    axes[2].set_ylabel("mismatch sum")
    axes[2].legend(loc="best")

    plt.tight_layout()
    out_path = out_dir / "nsbu23_route_capacity_hypothesis_dense.png"
    plt.savefig(out_path, dpi=150)
    print(f"saved {out_path}")

    data_out = out_dir / "nsbu23_route_capacity_hypothesis_dense_data.npz"
    meta_json = np.frombuffer(json.dumps(cfg, sort_keys=True).encode("utf-8"), dtype=np.uint8)
    savez_deterministic(
        data_out,
        meta_json=meta_json,
        start_idx=np.array(start_idx),
        r_steps=r_steps,
        p_best=np.array(-1 if p_best is None else p_best),
        max_violation=np.array(max_violation),
        worst_step=np.array(worst_step),
        e=e,
        e_tilde=e_tilde,
        partial_sum_e=partial_sum_e,
        psum_e_tilde=psum_e_tilde,
        c_tilde=c_tilde,
        verdict=np.array([verdict]),
        summability_status=np.array([summability_status]),
    )
    print(f"saved data: {data_out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
