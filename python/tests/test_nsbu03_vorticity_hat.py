import numpy as np


def test_vorticity_hat_sign() -> None:
    kx = np.array([1.0])
    ky = np.array([2.0])
    kz = np.array([3.0])
    u_hat = np.array([[4.0], [5.0], [6.0]], dtype=complex)

    omega_hat = np.zeros_like(u_hat)
    omega_hat[0] = 1j * (ky * u_hat[2] - kz * u_hat[1])
    omega_hat[1] = 1j * (kz * u_hat[0] - kx * u_hat[2])
    omega_hat[2] = 1j * (kx * u_hat[1] - ky * u_hat[0])

    expected = np.array([[-3j], [6j], [-3j]])
    assert np.allclose(omega_hat, expected)
