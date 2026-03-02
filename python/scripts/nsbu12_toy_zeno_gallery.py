from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "python" / "src"
sys.path.insert(0, str(SRC))

from nswave.toy_zeno import toy_case_series
from nswave.deterministic_npz import savez_deterministic


CASE_INFO = [
    ("Z1_fast_capacity", "Zeno due to fast-growing capacity => DIV fails"),
    ("Z2_vanishing_work", "Zeno due to vanishing work => WORK fails"),
    ("NZ_polynomial_capacity", "No-Zeno control => DIV holds (slow capacity growth)"),
]


def main() -> int:
    parser = argparse.ArgumentParser(description="NS-BU-12 toy Zeno gallery")
    parser.add_argument("--J", type=int, default=200, help="number of ladder levels")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out-dir", default=str(ROOT / "python" / "artifacts" / "paper"))
    args = parser.parse_args()

    J = args.J
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cfg = {
        "script": Path(__file__).name,
        "seed": args.seed,
        "N": 0,
        "dt": 0.0,
        "t_max": 0.0,
        "nu": 0.0,
        "mu": 0.0,
        "alpha": 0.0,
        "J": J,
        "out_dir": str(out_dir),
    }
    print(f"FIGURE_CONFIG: {json.dumps(cfg, sort_keys=True)}")
    j_plot = np.arange(1, J + 1)

    print(f"NS-BU-12 toy Zeno gallery (J={J})")
    series = {}
    for case_id, interp in CASE_INFO:
        data = toy_case_series(case_id, J)
        t_J = data["t_J"]
        sum_w_over_cap = data["sum_w_over_cap_J"]
        print(f"{case_id}: t_final={t_J[-1]:.6g} sum_w_over_cap={sum_w_over_cap[-1]:.6g}")
        print(f"interpretation: {interp}")
        series[case_id] = t_J

    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    for case_id, _ in CASE_INFO:
        ax.plot(j_plot, series[case_id], label=case_id)
    ax.set_xlabel("J")
    ax.set_ylabel("t_J (partial sum)")
    ax.set_title("Toy Zeno gallery: t_J vs J")
    ax.legend()
    ax.grid(True, alpha=0.3)

    out_path = out_dir / "nsbu12_toy_zeno_gallery.png"
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)

    print(f"saved {out_path.resolve()}")
    data_path = out_dir / "nsbu12_toy_zeno_gallery_data.npz"
    meta_json = np.frombuffer(json.dumps(cfg, sort_keys=True).encode("utf-8"), dtype=np.uint8)
    savez_deterministic(
        data_path,
        meta_json=meta_json,
        j_plot=j_plot,
        z1=series["Z1_fast_capacity"],
        z2=series["Z2_vanishing_work"],
        nz=series["NZ_polynomial_capacity"],
    )
    print(f"saved data: {data_path.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
