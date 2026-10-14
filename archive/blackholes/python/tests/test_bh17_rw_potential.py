import numpy as np

from bhwave.barrier import (
    r_of_rstar,
    regge_wheeler_V_rstar,
    rstar_of_r,
    sample_regge_wheeler_potential,
)


def test_rstar_inversion():
    for r in [2.5, 3.0, 5.0, 20.0]:
        x = rstar_of_r(r, M=1.0)
        r2 = r_of_rstar(x, M=1.0)
        assert abs(r2 - r) / r < 1e-10


def test_regge_wheeler_potential_tails():
    x = np.linspace(-40.0, 80.0, 2001)
    V = regge_wheeler_V_rstar(x, M=1.0, ell=2)
    assert float(np.max(V)) > 0.0
    assert float(V[0]) < 1e-3
    assert float(V[-1]) < 2e-3


def test_sample_helper_centers_peak():
    x, V, _ = sample_regge_wheeler_potential(shift_peak_to_zero=True)
    i_peak = int(np.argmax(V))
    assert abs(float(x[i_peak])) < 1e-9
