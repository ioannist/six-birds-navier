"""Analytic controls for the spectral solver and its original SGS port."""

import numpy as np
import pytest

from nswave.coarse_grain import restrict_hat_2d, restrict_hat_3d
from nswave.dyadic import dyadic_shell_mask
from nswave.ns2d_spectral import dealias_mask_2d, step_etdrk4
from nswave.ns3d_spectral import (
    _step_ifrk4,
    make_deterministic_u0_hat_3d,
    nonlinear_term_hat_3d,
)
from nswave.operators import dealias_mask_3d, kgrid_3d


@pytest.mark.parametrize("step_fn", [step_etdrk4, _step_ifrk4])
@pytest.mark.parametrize("quadratic", [False, True])
def test_integrating_factor_fourth_order_against_exact_solution(step_fn, quadratic):
    # y' = L*y + y^2, or y' = L*y + 1. Both have explicit solutions.
    linear = np.array([-0.7])
    y0, end = 0.6, 0.8
    exponential = np.exp(linear[0] * end)
    if quadratic:
        exact = y0 * exponential / (1 - y0 * (exponential - 1) / linear[0])
    else:
        exact = y0 * exponential + (exponential - 1) / linear[0]
    errors = []
    for nsteps in [4, 8, 16]:
        dt = end / nsteps
        y = np.array([y0])
        for n in range(nsteps):
            y = step_fn(y, dt, np.exp(linear * dt), np.exp(linear * dt / 2),
                        linear, lambda x, t: x * x if quadratic else np.ones_like(x),
                        n * dt)
        errors.append(abs(y[0] - exact))
    assert min(np.array(errors[:-1]) / np.array(errors[1:])) > 13


@pytest.mark.parametrize("N", [12, 15, 24])
def test_dealiasing_has_strict_quadratic_cutoff(N):
    modes = np.fft.fftfreq(N) * N
    for dim, mask_fn in [(2, dealias_mask_2d), (3, dealias_mask_3d)]:
        grids = np.meshgrid(*([modes] * dim), indexing="ij")
        mask = mask_fn(N).astype(bool)
        for grid in grids:
            assert 3 * np.max(np.abs(grid[mask])) < N


def test_dealiased_3d_nonlinearity_cancels_energy_at_divisible_by_three_grid():
    N = 12
    u_hat = make_deterministic_u0_hat_3d(N, seed=17, k_init=None)
    kx, ky, kz, k2 = kgrid_3d(N)
    nl = nonlinear_term_hat_3d(u_hat, kx, ky, kz, k2, mask=dealias_mask_3d(N))
    pairing = np.vdot(u_hat, nl).real
    assert abs(pairing) < 1e-13 * np.linalg.norm(u_hat) * np.linalg.norm(nl)


@pytest.mark.parametrize("dim", [2, 3])
def test_fourier_restriction_preserves_reality_and_nesting(dim):
    rng = np.random.default_rng(71)
    field = rng.standard_normal((16,) * dim)
    u_hat = np.fft.fftn(field)
    restrict = restrict_hat_2d if dim == 2 else restrict_hat_3d
    coarse = restrict(u_hat, 8)
    phys = np.fft.ifftn(coarse)
    assert np.max(np.abs(phys.imag)) < 1e-14
    assert np.allclose(restrict(u_hat, 4), restrict(coarse, 4), atol=1e-14)
    assert np.array_equal(restrict(u_hat, 16), u_hat)


def test_high_high_to_low_sgs_power_is_linear_in_small_shell_amplitude():
    # Exact trigonometric witness for the invalid shell-only bound in the old notes.
    N = 24
    x, y, _ = np.meshgrid(*([2 * np.pi * np.arange(N) / N] * 3), indexing="ij")
    a = np.array([0.0, 1.0, 0.0])[:, None, None, None]
    b = np.array([1.0, 3.0, 0.0])[:, None, None, None]
    c = np.array([1.0, -1.0, 0.0])[:, None, None, None]
    high = a * np.cos(4 * x) + b * np.cos(-3 * x + y)
    kx, ky, kz, k2 = kgrid_3d(N)
    shell = dyadic_shell_mask(np.sqrt(k2), 0)
    for eps in [0.1, 0.01, 0.001]:
        low = -eps * c * np.sin(x + y)
        full = high + low
        full_hat = np.fft.fftn(full, axes=(1, 2, 3))
        low_hat = np.fft.fftn(low, axes=(1, 2, 3))
        output_hat = (
            nonlinear_term_hat_3d(full_hat, kx, ky, kz, k2, mask=dealias_mask_3d(N))
            - nonlinear_term_hat_3d(low_hat, kx, ky, kz, k2, mask=dealias_mask_3d(N))
        ) * shell
        output = np.fft.ifftn(output_hat, axes=(1, 2, 3)).real
        assert np.allclose(output, -1.5 * c * np.sin(x + y), atol=1e-13)
        power = np.mean(np.sum(low * output, axis=0))
        shell_norm_sq = np.mean(np.sum(low * low, axis=0))
        assert np.isclose(power, 1.5 * eps, rtol=1e-12)
        assert np.isclose(power / shell_norm_sq, 1.5 / eps, rtol=1e-12)
