import numpy as np

from bhwave.kramers_kronig import kk_static_re0_from_im


def Z_debye_sum(omega: np.ndarray, a0: float, a_list: list[float], b_list: list[float]) -> np.ndarray:
    Z = np.full_like(omega, a0, dtype=np.complex128)
    for a_k, b_k in zip(a_list, b_list):
        Z = Z + a_k / (1.0 - 1j * omega / b_k)
    return Z


def test_static_kk_reconstruction():
    omega = np.logspace(-3, 3, 1201)
    a0 = 0.1
    a_list = [0.4, 0.3, 0.2]
    b_list = [0.8, 3.0, 15.0]

    Z = Z_debye_sum(omega, a0, a_list, b_list)
    ReZ0_true = a0 + sum(a_list)
    ReZ0_pred = kk_static_re0_from_im(omega, Z.imag, re_inf=a0)

    assert abs(ReZ0_pred - ReZ0_true) < 2e-2
