"""Controls for uncovered mathematical boundaries of the original experiments."""

import itertools
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from nswave.dyadic import percentile_frontier, shell_energy
from nswave.dyadic_models import simulate_kp
from nswave.coarse_grain import restrict_hat_3d
from nswave.galerkin_core import make_psd_operator
from nswave.illegal_initial_conditions import rescale_energy
from nswave.markov_lattice import init_sinusoidal_imbalance_ensemble, step_parallel_two_species_asep
from nswave.ns2d_spectral import simulate_ns2d
from nswave.ns3d_spectral import make_deterministic_u0_hat_3d, nonlinear_term_hat_3d, simulate_ns3d
from nswave.operators import dealias_mask_3d, kgrid_3d
from nswave.route_capacity_hypothesis import assess_mismatch_summability, classify_summability
from nswave.spt_legal_cert import compute_spt_legal_metrics


@pytest.mark.parametrize("horizon", [0.003, 0.1])
def test_dyadic_model_returns_the_requested_terminal_state(horizon):
    # J=1 has no cascade term: its independent exact solution is exponential.
    nu, force, initial = 0.7, 0.2, 1.0
    res = simulate_kp(J=1, dt=0.03, t_max=horizon, nu=nu, mu=0,
                      alpha=1.25, u0=np.array([initial]), forcing=np.array([force]),
                      record_every=20)
    expected = initial*np.exp(-nu*horizon) + force/nu*(-np.expm1(-nu*horizon))
    assert res["success"]
    assert res["t"][-1] == horizon
    assert res["u_hist"][-1, 0] == pytest.approx(expected, abs=1e-9)


def test_dyadic_nonfinite_diagnostic_cannot_report_success():
    with np.errstate(over="ignore", invalid="ignore"):
        res = simulate_kp(J=1, dt=0.01, t_max=0.1, nu=0, mu=0,
                          alpha=1.25, u0=np.array([1e200]))
    assert not res["success"]


@pytest.mark.parametrize("N", [2, 3, 4])
def test_scalar_parseval_uses_spatial_dimensions_even_on_small_grids(N):
    field = np.random.default_rng(17).standard_normal((N, N, N))
    got = shell_energy(np.fft.fftn(field), np.ones_like(field), L=2*np.pi)
    assert np.isclose(got, 0.5*np.mean(field**2)*(2*np.pi)**3)


def test_percentile_frontier_retains_the_shell_at_its_cutoff():
    assert percentile_frontier(np.array([1.0, 0, 0]), 0.01) == 0
    assert percentile_frontier(np.array([0.8, 0.15, 0.05]), 0.1) == 1
    assert percentile_frontier(np.array([0, 1.0, 0]), 0) == 1
    assert percentile_frontier(np.zeros(3), 0.1) == -1


def test_normalization_cannot_pay_raw_summability():
    n = np.arange(1, 41)
    result = classify_summability(np.ones(40), 1/n**2)
    assert result["status"] == "FAIL"
    assert result["normalized_status"] == "PASS"
    assert result["status_scope"] == "sampled_raw_trend_only"
    assert not assess_mismatch_summability(np.array([]), np.array([]), np.array([]))["status"].startswith("PASS")


def test_zero_capacity_prefix_has_no_step_factor_fit(tmp_path):
    data = tmp_path/"zero.npz"
    np.savez(data, C_tilde=np.zeros(4), e_route=np.zeros(4), e_route_tilde=np.zeros(4),
             partial_sum_e_route=np.zeros(4), U_rms=1.0, k_cut=np.ones(4), Lambda_direct=np.zeros(4))
    script = Path(__file__).resolve().parents[1]/"scripts/nsbu23_route_capacity_hypothesis_check.py"
    proc = subprocess.run([sys.executable, str(script), "--npz", str(data), "--out-dir", str(tmp_path)],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout+proc.stderr
    assert "OVERALL: INCONCLUSIVE" in proc.stdout
    assert "TAME_STEP_FACTOR: INCONCLUSIVE" in proc.stdout
    assert (tmp_path/"nsbu23_route_capacity_hypothesis_dense_data.npz").exists()


def test_disjoint_markov_matchings_conserve_each_species_exhaustively():
    # All 3^4 configurations, and all deterministic choices of the four probabilities.
    initial = np.array(list(itertools.product([-1, 0, 1], repeat=4)), dtype=np.int8)
    before = np.stack([np.sum(initial == species, axis=1) for species in [-1, 1]])
    for probabilities in itertools.product([0.0, 1.0], repeat=4):
        state = initial.copy()
        step_parallel_two_species_asep(state, *probabilities, np.random.default_rng(0))
        after = np.stack([np.sum(state == species, axis=1) for species in [-1, 1]])
        assert np.array_equal(before, after)
    with pytest.raises(ValueError, match="even lattice"):
        step_parallel_two_species_asep(np.array([[0, 1, 1]]), 1, 1, 0, 0, np.random.default_rng(0))
    with pytest.raises(ValueError, match="probabilities"):
        step_parallel_two_species_asep(initial.copy(), np.nan, 0, 0, 0, np.random.default_rng(0))
    with pytest.raises(ValueError, match="eps"):
        init_sinusoidal_imbalance_ensemble(4, 2, 1, 0.1, 0.2)


@pytest.mark.parametrize("dim", [2, 3])
@pytest.mark.parametrize("end", [0.07, 0.007])
def test_spectral_horizon_and_final_snapshot_agree_with_exact_shear(dim, end):
    N = 8
    x = 2*np.pi*np.arange(N)/N
    if dim == 2:
        omega = np.broadcast_to(np.sin(x)[:, None], (N, N))
        hat = np.fft.fftn(omega)
        result = simulate_ns2d(N=N, dt=0.03, t_max=end, nu=0.2, mu=0.1, alpha=1.25,
                              omega0_hat=hat, record_every=9, record_omega_hat_every=9)
        final, snapshots, t_snap = result.omega_hat_final, result.omega_hat_hist, result.t_omega
    else:
        velocity = np.zeros((3, N, N, N))
        velocity[1] = np.sin(x)[:, None, None]
        hat = np.fft.fftn(velocity, axes=(1, 2, 3))
        result = simulate_ns3d(N=N, dt=0.03, t_max=end, nu=0.2, mu=0.1, alpha=1.25,
                              u0_hat=hat, record_every=9, record_u_hat_every=9)
        final, snapshots, t_snap = result.u_hat_final, result.u_hat_hist, result.t_hat
    assert result.t[-1] == end and t_snap[-1] == end
    assert np.allclose(final, hat*np.exp(-0.3*end), atol=1e-13)
    assert np.array_equal(final, snapshots[-1])


def test_smooth_reference_recovery_rate_on_exact_nonlinearity_free_mode():
    N, end, nu, alpha = 8, 0.2, 0.1, 1.25
    x = 2*np.pi*np.arange(N)/N
    field = np.zeros((3, N, N, N))
    field[1] = np.sin(x)[:, None, None]
    initial = np.fft.fftn(field, axes=(1, 2, 3))
    reference = initial*np.exp(-nu*end)
    for mu in [1e-2, 1e-3, 1e-4]:
        res = simulate_ns3d(N=N, dt=0.03, t_max=end, nu=nu, mu=mu, alpha=alpha, u0_hat=initial)
        normalized_error = np.linalg.norm(res.u_hat_final-reference)/np.linalg.norm(reference)
        assert np.isclose(normalized_error, -np.expm1(-mu*end), rtol=1e-8)
        assert np.isclose(normalized_error/mu, end, rtol=0.002)


def test_nonreal_spectrum_cannot_be_dropped_into_a_passing_real_certificate():
    N = 8
    # Symmetric imaginary coefficients give a purely imaginary physical cosine.
    hist = np.zeros((1, 3, N, N, N), dtype=complex)
    hist[0, 1, 1, 0, 0] = hist[0, 1, -1, 0, 0] = 1j
    with pytest.raises(ValueError, match="Hermitian"):
        compute_spt_legal_metrics(hist, N=N, L=2*np.pi, window_c=1, energy_threshold=0.95)
    with pytest.raises(ValueError, match="Hermitian"):
        simulate_ns3d(N=N, dt=0.001, t_max=0.002, nu=0.1, mu=0, alpha=1.25, u0_hat=hist[0])


def test_aliased_control_still_preserves_a_real_solenoidal_spectral_state():
    N = 8
    initial = make_deterministic_u0_hat_3d(N, seed=29, k_init=None, dealias=False)
    res = simulate_ns3d(N=N, dt=0.001, t_max=0.002, nu=0.1, mu=0, alpha=1.25,
                       u0_hat=initial, dealias=False)
    physical = np.fft.ifftn(res.u_hat_final, axes=(1, 2, 3))
    assert np.max(np.abs(physical.imag)) < 1e-13
    assert np.max(res.div_rms) < 1e-13


def test_nonfinite_force_does_not_report_successful_simulation():
    N = 8
    with pytest.raises(ValueError, match="finite"):
        simulate_ns3d(N=N, dt=0.001, t_max=0.002, nu=0.1, mu=0, alpha=1.25,
                       forcing_hat_fn=lambda t, kx, ky, kz: np.full((3,N,N,N), np.nan))


def test_actual_nonlinear_two_stage_closure_retains_the_first_stage_response():
    hi, mid, co = 24, 12, 8
    u_hi = make_deterministic_u0_hat_3d(hi, seed=67, k_init=None)
    u_mid = restrict_hat_3d(u_hi, mid)*dealias_mask_3d(mid)
    u_co = restrict_hat_3d(u_hi, co)*dealias_mask_3d(co)
    nonlinearity = lambda u, n: nonlinear_term_hat_3d(u, *kgrid_3d(n), mask=dealias_mask_3d(n))
    nl_hi, nl_mid, nl_co = nonlinearity(u_hi, hi), nonlinearity(u_mid, mid), nonlinearity(u_co, co)
    direct = nl_co - restrict_hat_3d(nl_hi, co)
    second_only = nl_co - restrict_hat_3d(nl_mid, co)
    first_response = restrict_hat_3d(nl_mid - restrict_hat_3d(nl_hi, mid), co)
    assert np.linalg.norm(first_response) > 1e-6
    assert np.linalg.norm(direct-second_only) > 1e-6
    assert np.allclose(direct, second_only+first_response, atol=1e-12)


def test_ineligible_passive_fit_preserves_the_negative_observation(tmp_path):
    np.savez(tmp_path/"nsbu04_sgs_H_data.npz", omega=np.arange(6), H=-np.ones((1,6), dtype=complex),
             denom=np.ones((1,6)), shells=np.array([1]))
    script = Path(__file__).resolve().parents[1]/"scripts/nsbu05_sgs_realization_fit.py"
    proc = subprocess.run([sys.executable, str(script), "--out-dir", str(tmp_path)],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout+proc.stderr
    assert "NO SHELLS PASSED SELECTION" in proc.stdout
    with np.load(tmp_path/"nsbu05_sgs_fit_data.npz") as data:
        assert data["status"] == "NO_ELIGIBLE_SHELL"
        assert np.all(np.real(data["H"]) < 0)


def test_energy_matching_cannot_erase_a_tiny_nonzero_wavepacket():
    original = make_deterministic_u0_hat_3d(8, seed=73)
    for amplitude in [1.0, 1e-150]:
        matched = rescale_energy(amplitude*original, 2.0)
        physical = np.fft.ifftn(matched, axes=(1,2,3)).real
        energy = 0.5*np.mean(np.sum(physical**2, axis=0))*(2*np.pi)**3
        assert np.isclose(energy, 2.0, rtol=1e-12)
    with pytest.raises(ValueError, match="zero field"):
        rescale_energy(np.zeros_like(original), 1.0)


def test_psd_constructor_cannot_accept_an_antipassive_scale():
    with pytest.raises(ValueError, match="PSD scale"):
        make_psd_operator(3, scale=-1.0)


@pytest.mark.parametrize("script_name,data_name", [
    ("nsbu29_antilocalization_sweep.py", "nsbu29_antilocalization_sweep_data.npz"),
    ("ns_spt_legal_sweep.py", "ns_spt_legal_sweep_data.npz"),
])
def test_a_sweep_cannot_pass_by_ignoring_a_failed_configuration(tmp_path, script_name, data_name):
    script = Path(__file__).resolve().parents[1]/"scripts"/script_name
    proc = subprocess.run([sys.executable, str(script), "--resolutions", "8", "-1", "--nu-list", "0.001",
                           "--seeds", "0", "--dt", "0.001", "--t-max", "0.002",
                           "--record-u-hat-every", "1", "--out-dir", str(tmp_path)],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout+proc.stderr
    assert "SWEEP_STATUS: INCONCLUSIVE" in proc.stdout
    with np.load(tmp_path/data_name) as data:
        assert not data["coverage_complete"]
        assert data["failures"].size == 1
