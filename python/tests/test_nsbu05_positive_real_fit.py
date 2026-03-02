import numpy as np

from nswave.positive_real_fit import debye_model, fit_debye_positive_real


def test_positive_real_fit_synthetic() -> None:
    rng = np.random.default_rng(0)
    omega = np.logspace(-1, 1, 80)
    a0 = 0.2
    a = np.array([0.5, 0.3])
    b = np.array([0.8, 3.0])
    H_true = debye_model(omega, a0, a, b)
    noise = 0.01 * (rng.standard_normal(omega.shape) + 1j * rng.standard_normal(omega.shape))
    H_obs = H_true + noise

    res = fit_debye_positive_real(omega, H_obs, n_modes=2, max_nfev=300, seed=0)

    H_fit = debye_model(omega, res.a0, res.a, res.b)
    rmse = np.sqrt(np.mean(np.abs(H_fit - H_true) ** 2))
    assert res.success
    assert np.min(np.real(H_fit)) >= -1e-6
    assert rmse < 0.1 * np.median(np.abs(H_true))
