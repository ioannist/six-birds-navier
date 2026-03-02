import numpy as np

from nswave.capacity_estimator import capacity_L2_ratio, capacity_sliding_window


def test_capacity_exact_scalar_gain_positive() -> None:
    rng = np.random.default_rng(0)
    T, d = 2000, 5
    u = rng.standard_normal((T, d))
    lam = 1.7
    y = lam * u
    est = capacity_L2_ratio(u, y, dt=0.01)
    assert abs(est - lam) < 1e-10


def test_capacity_exact_scalar_gain_negative() -> None:
    rng = np.random.default_rng(1)
    T, d = 2000, 5
    u = rng.standard_normal((T, d))
    lam = -2.0
    y = lam * u
    est = capacity_L2_ratio(u, y, dt=0.01)
    assert abs(est) < 1e-12


def test_capacity_orthogonal_noise() -> None:
    rng = np.random.default_rng(2)
    T, d = 2000, 5
    u = rng.standard_normal((T, d))
    noise = rng.standard_normal((T, d))
    u_norm2 = np.sum(u * u, axis=1, keepdims=True)
    proj = np.sum(u * noise, axis=1, keepdims=True) / (u_norm2 + 1e-12)
    noise_orth = noise - proj * u
    lam = 0.9
    y = lam * u + 5.0 * noise_orth
    est = capacity_L2_ratio(u, y, dt=0.01)
    assert abs(est - lam) < 1e-10


def test_capacity_sliding_window_tracks_piecewise() -> None:
    rng = np.random.default_rng(3)
    T, d = 2000, 4
    u = rng.standard_normal((T, d))
    lam1, lam2 = 0.5, 2.0
    lam = np.ones(T) * lam1
    lam[T // 2 :] = lam2
    y = lam[:, None] * u
    t_mid, lam_series = capacity_sliding_window(u, y, dt=0.01, window_len=400, hop=200)
    assert t_mid.size == lam_series.size
    assert abs(lam_series[0] - lam1) < 1e-10
    assert abs(lam_series[-1] - lam2) < 1e-10


def test_capacity_zero_input_safe() -> None:
    u = np.zeros((100, 3))
    y = np.random.default_rng(4).standard_normal((100, 3))
    est = capacity_L2_ratio(u, y, dt=0.01)
    assert abs(est) < 1e-12
