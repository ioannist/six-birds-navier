import numpy as np

from bhwave.kramers_kronig import kk_reconstruct_re_from_im


def debye_sum(omega: np.ndarray, a0: float, a_list: list[float], b_list: list[float]) -> np.ndarray:
    Z = np.full_like(omega, a0, dtype=np.complex128)
    for a_k, b_k in zip(a_list, b_list):
        Z = Z + a_k / (1.0 - 1j * omega / b_k)
    return Z


def test_kk_reconstruction_case1():
    omega_grid = np.logspace(-2, 2, 301)
    a0 = 0.2
    a_list = [0.7]
    b_list = [3.0]

    Z_true = debye_sum(omega_grid, a0, a_list, b_list)
    re_pred = kk_reconstruct_re_from_im(omega_grid, Z_true.imag, re_inf=a0)

    n = omega_grid.size
    idx0 = int(0.2 * n)
    idx1 = int(0.8 * n)
    max_abs_err = np.max(np.abs(re_pred[idx0:idx1] - Z_true.real[idx0:idx1]))

    assert max_abs_err < 5e-2
