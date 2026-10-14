import numpy as np

from nswave.ns2d_spectral import make_deterministic_omega0_hat


def test_real_initial_condition_ifft2():
    N = 32
    omega_hat = make_deterministic_omega0_hat(N, seed=0, k_init=10, dealias=True)
    omega = np.fft.ifft2(omega_hat)

    max_im = float(np.max(np.abs(omega.imag)))
    max_re = float(np.max(np.abs(omega.real)))
    assert max_im <= 1e-12 * max(1.0, max_re)
