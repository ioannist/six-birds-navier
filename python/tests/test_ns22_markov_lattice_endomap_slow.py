import numpy as np
import pytest

from nswave.markov_lattice import (
    init_sinusoidal_imbalance_ensemble,
    step_parallel_two_species_asep,
)


@pytest.mark.slow
def test_ns22_markov_endomap_scaling_slow() -> None:
    N = 64
    M = 1000
    n_steps = 300
    rho0 = 0.05
    eps = 0.02
    k_list = [1, 2, 3]

    def run_case(
        p_plus: float,
        q_plus: float,
        p_minus: float,
        q_minus: float,
        k_fit_count: int,
        seeds: list[int],
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
        k_phys2 = []
        lam_list = []
        v_list = []
        snapshot_steps = np.linspace(0, n_steps, 9, dtype=int).tolist()
        snapshot_steps = sorted(set(snapshot_steps))

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
                            state, p_plus, q_plus, p_minus, q_minus, rng
                        )

                amps_seed = np.array(amps, dtype=complex)
                times_arr = np.array(times, dtype=float)
                amps_stack.append(amps_seed)
                phase_seed = np.unwrap(np.angle(amps_seed))
                phase_slope_seed, _ = np.polyfit(times_arr, phase_seed, 1)
                v_vals.append(phase_slope_seed / k_phys if k_phys > 0 else 0.0)

            amps_stack = np.array(amps_stack)
            amps_mean = np.mean(amps_stack, axis=0)
            times = times_arr
            amp_abs = np.abs(amps_mean)
            amp0 = amp_abs[0]
            valid = amp_abs >= 0.2 * amp0
            if np.count_nonzero(valid) < 3:
                valid = np.arange(len(times)) < 3

            slope, _ = np.polyfit(times[valid], np.log(amp_abs[valid]), 1)
            lam = -float(slope)
            D_est = lam / (k_phys**2)
            v_est = float(np.median(np.abs(v_vals)))

            k_phys2.append(k_phys**2)
            lam_list.append(lam)
            v_list.append(v_est)
            assert lam > 0
            assert D_est > 0

        k_phys2 = np.array(k_phys2)
        lam_arr = np.array(lam_list)
        v_arr = np.array(v_list)
        k_fit = k_phys2[:k_fit_count]
        lam_fit = lam_arr[:k_fit_count]
        D_fit = float(np.sum(k_fit * lam_fit) / np.sum(k_fit**2))
        lam_pred = D_fit * k_fit
        ss_res = float(np.sum((lam_fit - lam_pred) ** 2))
        ss_tot = float(np.sum((lam_fit - lam_fit.mean()) ** 2))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
        return k_phys2, lam_arr, v_arr, r2

    seeds = [0, 1, 2]
    k_phys2, lam_eq, v_eq, r2_eq = run_case(0.5, 0.5, 0.5, 0.5, 3, seeds)
    _, lam_dr, v_dr, r2_dr = run_case(0.54, 0.46, 0.46, 0.54, 3, seeds)

    assert r2_eq >= 0.75
    assert r2_dr >= 0.70
    assert abs(v_dr[0]) >= 1e-2
    assert np.all(lam_eq > 0)
    assert np.all(lam_dr > 0)
