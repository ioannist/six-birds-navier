import os
import sys
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "python" / "artifacts"
LOGS = ART / "logs"
SCRIPTS_WITH_PAPER_ARGS = {
    "ns14_energy_ledger_3d.py",
    "ns15_uv_completion_turns_on_3d.py",
    "ns16_mu_to_zero_convergence_3d.py",
    "ns18_exponent_sanity_harness.py",
}

SCRIPTS = [
    "ns05_energy_ledger.py",
    "ns06_uv_completion_turns_on.py",
    "ns07_mu_to_zero_convergence.py",
    "ns08_infer_closure_kernel.py",
    "ns09_invariance_sanity.py",
    "ns13_smoke_demo.py",
    "ns14_energy_ledger_3d.py",
    "ns15_uv_completion_turns_on_3d.py",
    "ns16_mu_to_zero_convergence_3d.py",
    "ns18_exponent_sanity_harness.py",
]

ARTIFACTS = [
    "ns05_energy_ledger.png",
    "ns06_uv_completion.png",
    "ns07_mu_convergence.png",
    "ns08_closure_kernel.png",
    "ns09_invariance.png",
    "ns14_energy_ledger_3d.png",
    "ns15_uv_completion_3d.png",
    "ns16_mu_convergence_3d.png",
    "ns18_exponent_sanity.png",
]

INCLUDE_SLOW = os.environ.get("NS_INCLUDE_SLOW", "") == "1"
if INCLUDE_SLOW:
    SCRIPTS.append("ns20_infer_closure_kernel_3d.py")
    ARTIFACTS.append("ns20_closure_kernel_3d.png")
    SCRIPTS.append("ns21_infer_closure_kernel_3d_forced_avg.py")
    ARTIFACTS.append("ns21_closure_kernel_3d_forced_avg.png")

INCLUDE_OPTIONAL = os.environ.get("NS_INCLUDE_OPTIONAL", "") == "1"
if INCLUDE_OPTIONAL:
    SCRIPTS.append("ns22_markov_lens_endomap_demo.py")
    ARTIFACTS.append("ns22_markov_laplacian_scaling.png")


def run_script(script: str) -> tuple[int, str]:
    script_path = ROOT / "python" / "scripts" / script
    cmd = [sys.executable, str(script_path)]
    if script in SCRIPTS_WITH_PAPER_ARGS:
        cmd.extend(["--seed", "0", "--out-dir", str(ART)])
    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
    )
    output = proc.stdout + proc.stderr
    return proc.returncode, output


def main() -> int:
    ART.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)

    results: dict[str, str] = {}
    ok = True

    for script in SCRIPTS:
        code, output = run_script(script)
        log_path = LOGS / f"{Path(script).stem}.log"
        log_path.write_text(output)
        status = "OK" if code == 0 else "FAIL"
        results[script] = status
        if code != 0:
            ok = False

    missing = []
    for art in ARTIFACTS:
        if not (ART / art).exists():
            missing.append(art)

    if missing:
        ok = False

    print("NS-12 run-all summary:")
    for script in SCRIPTS:
        print(f"  {script}: {results[script]}")
    if not INCLUDE_SLOW:
        print("  ns20_infer_closure_kernel_3d.py: SKIP (slow; set NS_INCLUDE_SLOW=1)")
        print("  ns21_infer_closure_kernel_3d_forced_avg.py: SKIP (slow; set NS_INCLUDE_SLOW=1)")
    if not INCLUDE_OPTIONAL:
        print("  ns22_markov_lens_endomap_demo.py: SKIP (optional; set NS_INCLUDE_OPTIONAL=1)")
    print(f"Logs directory: {LOGS}")

    if missing:
        print("Missing artifacts:")
        for art in missing:
            print(f"  {art}")
    else:
        print("Artifacts found:")
        for art in ARTIFACTS:
            print(f"  {art}")

    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
