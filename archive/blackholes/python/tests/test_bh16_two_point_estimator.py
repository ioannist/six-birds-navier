import numpy as np

from bhwave.inference import estimate_Rout_from_two_point_spectra


def test_two_point_estimator_recovers_constant_R():
    N = 512
    dt = 0.25
    omega_pos = 2.0 * np.pi * np.fft.rfftfreq(N, dt)

    y1 = 28.0
    y2 = 28.02

    R_true = 0.3 * np.exp(0.2j)
    A_in = np.exp(-0.5 * ((omega_pos - 1.5) / 0.5) ** 2)
    if A_in.size > 0:
        A_in[0] = 0.0
    A_out = R_true * A_in

    Psi1 = A_in * np.exp(-1j * omega_pos * y1) + A_out * np.exp(1j * omega_pos * y1)
    Psi2 = A_in * np.exp(-1j * omega_pos * y2) + A_out * np.exp(1j * omega_pos * y2)

    psi1 = np.fft.irfft(Psi1, n=N)
    psi2 = np.fft.irfft(Psi2, n=N)

    Psi1p = np.fft.rfft(psi1)
    Psi2p = np.fft.rfft(psi2)

    R_est, mask = estimate_Rout_from_two_point_spectra(
        omega_pos=omega_pos,
        Psi1=Psi1p,
        Psi2=Psi2p,
        y1=y1,
        y2=y2,
    )

    diff = np.abs(R_est[mask] - R_true)
    assert diff.size > 0
    assert float(np.max(diff)) < 1e-6
