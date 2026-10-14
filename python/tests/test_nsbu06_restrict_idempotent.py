import numpy as np

from nswave.coarse_grain import restrict_hat_3d


def test_restrict_idempotent() -> None:
    N_hi = 16
    N_mid = 8
    N_co = 4
    rng = np.random.default_rng(0)
    u_hat = rng.standard_normal((3, N_hi, N_hi, N_hi)) + 1j * rng.standard_normal(
        (3, N_hi, N_hi, N_hi)
    )

    u_co_1 = restrict_hat_3d(u_hat, N_co)
    u_co_2 = restrict_hat_3d(restrict_hat_3d(u_hat, N_mid), N_co)

    denom = max(1.0, np.max(np.abs(u_co_1)))
    assert np.max(np.abs(u_co_1 - u_co_2)) <= 1e-12 * denom
