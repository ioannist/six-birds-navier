#!/usr/bin/env python3
from __future__ import annotations

import pathlib
import subprocess
import sys


def run_script(script: pathlib.Path, logs_dir: pathlib.Path) -> tuple[bool, pathlib.Path]:
    log_path = logs_dir / f"{script.stem}.log"
    result = subprocess.run(
        [sys.executable, str(script)],
        capture_output=True,
        text=True,
    )
    log_path.write_text(result.stdout + result.stderr, encoding="utf-8")
    return result.returncode == 0, log_path


def main() -> int:
    root = pathlib.Path(__file__).resolve().parents[2]
    artifacts_dir = root / "python" / "artifacts"
    logs_dir = artifacts_dir / "logs"
    artifacts_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    scripts = [
        root / "python" / "scripts" / "bh02_convention_sanity.py",
        root / "python" / "scripts" / "bh03_free_case_demo.py",
        root / "python" / "scripts" / "bh04_passivity_checks.py",
        root / "python" / "scripts" / "bh05_kramers_kronig_numeric.py",
        root / "python" / "scripts" / "bh06_barrier_table.py",
        root / "python" / "scripts" / "bh17_rw_barrier_compare.py",
        root / "python" / "scripts" / "bh18_constraint_stacking.py",
        root / "python" / "scripts" / "bh07_echo_formula_validation.py",
        root / "python" / "scripts" / "bh15_echo_freqdep_validation.py",
        root / "python" / "scripts" / "bh08_time_domain_demo.py",
        root / "python" / "scripts" / "bh09_passive_fit_demo.py",
        root / "python" / "scripts" / "bh16_fdtd_inference_demo.py",
        root / "python" / "scripts" / "bh10_cross_channel_linkage.py",
        root / "python" / "scripts" / "bh11_superradiance_ledger_demo.py",
        root / "python" / "scripts" / "bh13_energy_ledger_demo.py",
    ]

    results: dict[str, bool] = {}
    logs: dict[str, pathlib.Path] = {}
    for script in scripts:
        ok, log_path = run_script(script, logs_dir)
        results[script.name] = ok
        logs[script.name] = log_path

    print("BH-12 run-all summary:")
    for name in scripts:
        status = "OK" if results[name.name] else "FAIL"
        print(f"  {name.name}: {status}")
    print(f"Logs directory: {logs_dir}")

    expected = [
        "bh05_kk_case2.png",
        "bh08_waveform.png",
        "bh09_fit_waveform.png",
        "bh17_rw_vs_pt.png",
        "bh10_cross_channel.png",
        "bh11_superradiance_gain.png",
        "bh13_energy.png",
        "bh15_freqdep_echo_error.png",
        "bh16_fdtd_fit_waveform.png",
        "bh16_fdtd_fit_R_error.png",
        "bh18_constraint_stacking.png",
    ]
    missing = [name for name in expected if not (artifacts_dir / name).exists()]
    if missing:
        print("Missing artifacts:")
        for name in missing:
            print(f"  {name}")
        return 1

    if not all(results.values()):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
