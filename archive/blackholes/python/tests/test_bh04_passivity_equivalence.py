import numpy as np

from bhwave.conventions import R_from_Z


def flux(omega: float, R: complex) -> float:
    return omega * (abs(R) ** 2 - 1.0)


def test_grid_sanity():
    for Z in [1 + 0j, 0.2 + 0.3j, 0 + 0.5j]:
        R = R_from_Z(Z)
        assert abs(R) <= 1 + 1e-12

    Z = -0.2 + 0j
    R = R_from_Z(Z)
    assert abs(R) >= 1 + 1e-3


def test_random_passive_samples():
    rng = np.random.default_rng(0)
    re = rng.uniform(0.05, 3.0, 200)
    im = rng.uniform(-3.0, 3.0, 200)
    Z_vals = re + 1j * im

    for Z in Z_vals:
        R = R_from_Z(Z)
        assert abs(R) <= 1 + 1e-12


def test_random_active_samples():
    rng = np.random.default_rng(0)
    re = rng.uniform(-3.0, -0.2, 200)
    im = rng.uniform(-3.0, 3.0, 200)
    Z_vals = re + 1j * im

    for Z in Z_vals:
        R = R_from_Z(Z)
        assert abs(R) >= 1 - 1e-12
        assert abs(R) > 1 + 1e-6


def test_flux_sign_check():
    omega = 2.0
    R_passive = R_from_Z(0.2 + 0.3j)
    R_active = R_from_Z(-0.2 + 0j)

    assert flux(omega, R_passive) <= 1e-12
    assert flux(omega, R_active) > 1e-3
