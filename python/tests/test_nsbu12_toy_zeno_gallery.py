import numpy as np

from nswave.toy_zeno import toy_case_series


def test_toy_zeno_gallery_bounds_and_monotone() -> None:
    J = 200
    cases = [
        "Z1_fast_capacity",
        "Z2_vanishing_work",
        "NZ_polynomial_capacity",
    ]

    for case_id in cases:
        data = toy_case_series(case_id, J)
        t_J = data["t_J"]
        assert t_J.shape == (J,)
        assert np.all(np.diff(t_J) >= -1e-12)

    t_z1 = toy_case_series("Z1_fast_capacity", J)["t_J"][-1]
    t_z2 = toy_case_series("Z2_vanishing_work", J)["t_J"][-1]
    t_nz = toy_case_series("NZ_polynomial_capacity", J)["t_J"][-1]

    assert abs(t_z1 - 2.0) < 1e-3
    assert abs(t_z2 - 2.0) < 1e-3
    assert t_nz > 5.0
