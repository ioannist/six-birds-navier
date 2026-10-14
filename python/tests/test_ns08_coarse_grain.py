import numpy as np

from nswave.coarse_grain import restrict_hat_2d


def test_restrict_hat_preserves_mode_amplitude():
    L = 2 * np.pi
    N_hi = 32
    N_co = 16

    kx = 3
    ky = 2

    omega_hat_hi = np.zeros((N_hi, N_hi), dtype=complex)
    amp = (N_hi**2) / 2.0

    ix = kx % N_hi
    iy = ky % N_hi
    ixm = (-kx) % N_hi
    iym = (-ky) % N_hi

    omega_hat_hi[ix, iy] = amp
    omega_hat_hi[ixm, iym] = amp

    omega_hat_co = restrict_hat_2d(omega_hat_hi, N_co)
    omega_co = np.fft.ifft2(omega_hat_co).real

    x = L * np.arange(N_co) / N_co
    y = L * np.arange(N_co) / N_co
    xx, yy = np.meshgrid(x, y, indexing="ij")
    omega_exact = np.cos(kx * xx + ky * yy)

    err = np.max(np.abs(omega_co - omega_exact))
    assert err < 1e-10
