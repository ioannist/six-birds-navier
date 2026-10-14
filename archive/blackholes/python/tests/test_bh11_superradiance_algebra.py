import cmath
import numpy as np

from bhwave.superradiance import growth_rate_from_gain, pole_estimate_from_G


def test_pole_estimate_from_gain():
    Dx = 8.0
    r = 0.96 * cmath.exp(0.1j)
    R0 = 1.2 * cmath.exp(-0.2j)
    G = R0 * r

    pole = pole_estimate_from_G(G, Dx, n=1)
    assert pole.imag > 0

    expected_im = np.log(abs(G)) / (2.0 * Dx)
    assert abs(pole.imag - expected_im) < 1e-12
    assert abs(growth_rate_from_gain(abs(G), Dx) - expected_im) < 1e-12
