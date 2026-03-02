import numpy as np

from bhwave.inference import (
    R_pred_grid,
    Z_debye_k1,
    estimate_Rpred_from_baseline_ratio,
    gaussian_source,
)


def test_ratio_estimator_recovers_R_true():
    N = 512
    dt = 0.25
    omega_pos = 2.0 * np.pi * np.fft.rfftfreq(N, dt)

    r = 0.4
    t = np.sqrt(1.0 - r * r)
    Dx = 8.0

    R_base = r * np.exp(-2j * omega_pos * Dx)

    Z_true = Z_debye_k1(omega_pos, 0.2, 0.5, 2.0)
    R0_true = (1.0 - Z_true) / (1.0 + Z_true)
    R_true = R_pred_grid(
        omega_pos,
        R0_true,
        r * np.ones_like(omega_pos, dtype=complex),
        t * np.ones_like(omega_pos, dtype=complex),
        Dx,
    )

    S_pos = gaussian_source(omega_pos, omega0=1.5, sigma=0.5)
    S_pos = S_pos.astype(complex)
    if S_pos.size > 0:
        S_pos[0] = 0.0

    F_base = S_pos * R_base
    F_true = S_pos * R_true
    y_base = np.fft.irfft(F_base, n=N)
    y_true = np.fft.irfft(F_true, n=N)

    R_est, mask = estimate_Rpred_from_baseline_ratio(
        waveform_true=y_true,
        waveform_base=y_base,
        R_base=R_base,
    )

    diff = np.abs(R_est[mask] - R_true[mask])
    assert diff.size > 0
    assert float(np.max(diff)) < 1e-8
