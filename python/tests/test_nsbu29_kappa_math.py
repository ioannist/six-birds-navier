import numpy as np

from nswave.localization import localization_factors, window_for_shell


def test_kappa_is_one_for_constructed_profile():
    N = 32
    j_values = np.array([1, 2, 3, 4], dtype=int)
    windows = np.array([window_for_shell(N, int(j), c=1.0) for j in j_values], dtype=int)
    baseline = (windows.astype(float) ** 3) / float(N**3)

    mean_conc = (j_values.astype(float) + 1.0) * baseline
    _, _, kappa = localization_factors(
        mean_conc,
        N=N,
        windows=windows,
        j_values=j_values,
    )

    assert np.allclose(kappa, np.ones_like(kappa), atol=1e-12, rtol=0.0)
