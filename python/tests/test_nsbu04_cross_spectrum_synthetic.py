import numpy as np


def test_cross_spectrum_single_mode() -> None:
    T = 128
    dt = 0.1
    freq0 = 1.0 / (T * dt)
    omega0 = 2 * np.pi * freq0
    lam = 0.2 + 0.1j

    t = np.arange(T) * dt
    u = np.exp(1j * omega0 * t)
    c = lam * u

    u_series = u[:, None]
    c_series = c[:, None]
    U = np.fft.fft(u_series, axis=0)
    C = np.fft.fft(c_series, axis=0)
    numer = C[:, 0] * np.conj(U[:, 0])
    denom = np.abs(U[:, 0]) ** 2
    H = numer / denom

    assert np.isclose(H[1], lam, atol=1e-2)
    assert np.isclose(np.real(H[1]), np.real(lam), atol=1e-2)
