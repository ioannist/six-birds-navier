import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "python" / "src"
sys.path.insert(0, str(SRC))

from nswave.markov_lattice import (  # noqa: E402
    init_sinusoidal_imbalance_ensemble,
    step_parallel_two_species_asep,
)


ART = ROOT / "python" / "artifacts"


def main() -> int:
    N = 64
    M = 1000
    n_steps = 300
    rho0 = 0.05
    eps = 0.02
    p_plus = 0.50
    q_plus = 0.50
    p_minus = 0.50
    q_minus = 0.50
    p_plus_drive = 0.54
    q_plus_drive = 0.46
    p_minus_drive = 0.46
    q_minus_drive = 0.54
    k_list = [1, 2, 3, 4, 5, 6]
    seeds = [0, 1, 2]

    ART.mkdir(parents=True, exist_ok=True)

    print("NS-22 Markov lens endomap (1D ASEP)")
    print(f"N={N} M={M} steps={n_steps} rho0={rho0} eps={eps}")

    snapshot_steps = np.linspace(0, n_steps, 9, dtype=int).tolist()
    snapshot_steps = sorted(set(snapshot_steps))

    def run_case(
        label: str,
        p_p: float,
        q_p: float,
        p_m: float,
        q_m: float,
        k_fit_count: int,
    ) -> tuple[np.ndarray, np.ndarray, float, float]:
        k_phys2 = []
        lam_list = []
        v_list = []
        print(label)
        print(f"p_plus={p_p} q_plus={q_p} p_minus={p_m} q_minus={q_m}")
        print(f"a_plus={np.log(p_p/q_p):.3f} a_minus={np.log(p_m/q_m):.3f}")
        print("lens: momentum (plus/minus imbalance)")
        print("k  k_phys^2   lam       D_est     v_est")

        for k_mode in k_list:
            amps_stack = []
            v_vals = []
            times_arr = None
            k_phys = 2 * np.pi * k_mode / N
            for seed in seeds:
                state = init_sinusoidal_imbalance_ensemble(
                    N, M, k_mode, rho0, eps, seed=seed
                )
                rng = np.random.default_rng(seed)
                amps = []
                times = []

                for step in range(n_steps + 1):
                    if step in snapshot_steps:
                        m = state.mean(axis=0)
                        mhat = np.fft.fft(m)
                        amps.append(mhat[k_mode])
                        times.append(step)
                    if step < n_steps:
                        step_parallel_two_species_asep(
                            state, p_p, q_p, p_m, q_m, rng
                        )

                amps_seed = np.array(amps, dtype=complex)
                times_arr = np.array(times, dtype=float)
                amps_stack.append(amps_seed)
                phase_seed = np.unwrap(np.angle(amps_seed))
                phase_slope_seed, _ = np.polyfit(
                    times_arr, phase_seed, 1
                )
                v_vals.append(phase_slope_seed / k_phys if k_phys > 0 else 0.0)

            amps_stack = np.array(amps_stack)
            amps_mean = np.mean(amps_stack, axis=0)
            times = times_arr
            amp_abs = np.abs(amps_mean)

            k_phys = 2 * np.pi * k_mode / N
            amp0 = amp_abs[0]
            valid = amp_abs >= 0.2 * amp0
            if np.count_nonzero(valid) < 3:
                valid = np.arange(len(times)) < 3
            slope, _ = np.polyfit(times[valid], np.log(amp_abs[valid]), 1)
            lam = -float(slope)
            phase = np.unwrap(np.angle(amps_mean))
            phase_slope, _ = np.polyfit(times[valid], phase[valid], 1)
            v_est = float(np.median(np.abs(v_vals)))

            D_est = lam / (k_phys**2)

            k_phys2.append(k_phys**2)
            lam_list.append(lam)
            v_list.append(v_est)

            print(
                f"{k_mode:>1d}  {k_phys**2:9.6f}  {lam:8.6f}  {D_est:8.6f}  {v_est:8.6f}"
            )

        k_phys2 = np.array(k_phys2)
        lam_arr = np.array(lam_list)

        k_fit = k_phys2[:k_fit_count]
        lam_fit = lam_arr[:k_fit_count]
        D_fit = float(np.sum(k_fit * lam_fit) / np.sum(k_fit**2))
        lam_pred = D_fit * k_fit
        ss_res = float(np.sum((lam_fit - lam_pred) ** 2))
        ss_tot = float(np.sum((lam_fit - lam_fit.mean()) ** 2))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

        print(f"D_fit={D_fit:.6f} R2={r2:.3f}")
        print(
            f"Note: D_fit/R2 computed on low-k subset (k=1..{k_fit_count}); fit_k_max={k_fit_count}; high-k deviations expected for lattice dynamics."
        )
        return k_phys2, lam_arr, D_fit, r2

    k_phys2_eq, lam_eq, D_fit_eq, r2_eq = run_case(
        "Equilibrium diffusion (p=q)", p_plus, q_plus, p_minus, q_minus, 3
    )
    k_phys2_dr, lam_dr, D_fit_dr, r2_dr = run_case(
        "Driven case (nonzero affinity)",
        p_plus_drive,
        q_plus_drive,
        p_minus_drive,
        q_minus_drive,
        3,
    )

    plt.figure(figsize=(6, 4))
    plt.plot(k_phys2_eq, lam_eq, "o", label="eq data")
    plt.plot(k_phys2_dr, lam_dr, "s", label="drive data")
    k_line = np.linspace(0, k_phys2_eq[3], 50)
    plt.plot(k_line, D_fit_eq * k_line, "--", label="eq fit")
    plt.plot(k_line, D_fit_dr * k_line, "--", label="drive fit")
    plt.xlabel(r"$k_{phys}^2$")
    plt.ylabel(r"$\lambda(k)$")
    plt.title("NS-22 Markov lens: Laplacian scaling (eq + drive)")
    plt.legend()
    plt.tight_layout()

    out_path = ART / "ns22_markov_laplacian_scaling.png"
    plt.savefig(out_path, dpi=150)
    print(f"saved {out_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
