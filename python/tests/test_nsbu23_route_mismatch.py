import numpy as np

from nswave.capacity_estimator import capacity_L2_ratio


def test_e_route_zero_for_telescoping():
    rng = np.random.default_rng(0)
    u = rng.standard_normal((100, 4))
    y_sum = rng.standard_normal((100, 4))
    y_direct = y_sum.copy()
    y_route = y_direct - y_sum
    lam = capacity_L2_ratio(u, y_route, dt=0.1)
    assert abs(lam) <= 1e-12


def test_e_route_positive_for_mismatch():
    rng = np.random.default_rng(1)
    u = rng.standard_normal((100, 4))
    y_sum = rng.standard_normal((100, 4))
    delta = rng.standard_normal((100, 4)) * 0.1
    y_direct = y_sum + delta
    y_route = y_direct - y_sum
    lam = capacity_L2_ratio(u, y_route, dt=0.1)
    assert lam > 0.0
