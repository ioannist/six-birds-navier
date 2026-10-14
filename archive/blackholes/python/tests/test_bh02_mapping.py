import math
import random

import pytest

from bhwave.conventions import (
    R_from_Y,
    R_from_Z,
    Y_from_R,
    Z_from_R,
    is_passive_R,
    is_passive_Z,
)


def test_sanity_cases():
    assert R_from_Z(1 + 0j) == pytest.approx(0 + 0j)
    assert R_from_Z(0 + 0j) == pytest.approx(1 + 0j)

    R = R_from_Z(0 + 0.5j)
    assert abs(abs(R) - 1.0) <= 1e-12

    Z = 0.2 + 0.3j
    R = R_from_Z(Z)
    assert abs(R) < 1.0
    assert is_passive_Z(Z)
    assert is_passive_R(R)

    Z = -0.2 + 0j
    R = R_from_Z(Z)
    assert abs(R) > 1.0
    assert not is_passive_Z(Z)
    assert not is_passive_R(R)


def test_round_trip_random_values():
    rng = random.Random(0)
    samples = []
    while len(samples) < 10:
        re = 0.1 + 2.0 * rng.random()
        im = -1.0 + 2.0 * rng.random()
        Z = re + 1j * im
        if abs(Z + 1) <= 0.1:
            continue
        samples.append(Z)

    for Z in samples:
        Z2 = Z_from_R(R_from_Z(Z))
        assert abs(Z2 - Z) <= 1e-12


def test_y_r_mapping_checks():
    omega = 2.0
    assert Y_from_R(0 + 0j, omega) == pytest.approx(-1j * omega)
    assert Y_from_R(1 + 0j, omega) == pytest.approx(0 + 0j)
    assert R_from_Y(-1j * omega, omega) == pytest.approx(0 + 0j)
    assert R_from_Y(0 + 0j, omega) == pytest.approx(1 + 0j)


def test_r_round_trip_random_values():
    rng = random.Random(0)
    omega = 2.0
    for _ in range(50):
        radius = 0.95 * rng.random()
        angle = 2.0 * math.pi * rng.random()
        R = radius * (math.cos(angle) + 1j * math.sin(angle))
        R2 = R_from_Y(Y_from_R(R, omega), omega)
        assert abs(R2 - R) <= 1e-12
