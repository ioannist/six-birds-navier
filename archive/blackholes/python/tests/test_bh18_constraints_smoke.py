import numpy as np


def test_constraints_smoke_algebra():
    a0 = 0.2
    a2 = 0.3
    b2 = 5.0

    a1_grid = np.linspace(0.0, 1.2, 5)
    b1_grid = np.logspace(np.log10(0.2), np.log10(5.0), 5)
    omega = np.linspace(0.2, 2.0, 20)

    a1 = a1_grid[:, None, None]
    b1 = b1_grid[None, :, None]
    om = omega[None, None, :]

    Z = a0 + a1 / (1.0 - 1j * om / b1) + a2 / (1.0 - 1j * om / b2)
    R0 = (1.0 - Z) / (1.0 + Z)

    assert float(np.min(np.real(Z))) >= -1e-12
    assert float(np.max(np.abs(R0))) <= 1.0 + 1e-12

    L_proxy = np.abs((a0 + a1_grid + a2) - 1.0)
    expected = np.abs((a0 + a1_grid + a2) - 1.0)
    assert np.allclose(L_proxy, expected)
