import numpy as np

from nswave.localization import max_local_energy_fraction


def test_max_local_energy_fraction_uniform_field():
    N = 16
    window = 4
    u = np.ones((3, N, N, N), dtype=float)
    frac = max_local_energy_fraction(u, window)

    expected = (window**3) / float(N**3)
    assert abs(frac - expected) <= 1e-10


def test_max_local_energy_fraction_spike_field():
    N = 16
    window = 2
    u = np.zeros((3, N, N, N), dtype=float)
    u[0, 0, 0, 0] = 1.0

    frac = max_local_energy_fraction(u, window)
    assert frac > 0.9
