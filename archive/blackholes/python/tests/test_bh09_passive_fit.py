import numpy as np

from bhwave.inference import (
    R_pred_grid,
    Z_debye_k1,
    fit_debye_k1_from_Rpred,
    passivity_residuals,
)


def test_fit_debye_k1_simple():
    omega_pos = np.linspace(0.2, 5.0, 80)
    r = 0.4 * np.ones_like(omega_pos, dtype=complex)
    t = np.sqrt(1.0 - r * r) * np.ones_like(omega_pos, dtype=complex)
    Dx = 8.0

    a0_true = 0.2
    a1_true = 0.5
    b1_true = 2.0

    Z_true = Z_debye_k1(omega_pos, a0_true, a1_true, b1_true)
    R0_true = (1.0 - Z_true) / (1.0 + Z_true)
    Rpred_target = R_pred_grid(omega_pos, R0_true, r, t, Dx)

    fit = fit_debye_k1_from_Rpred(
        omega_pos=omega_pos,
        Rpred_target=Rpred_target,
        r=r,
        t=t,
        Dx=Dx,
        p0=(0.15, 0.3, 1.0),
    )

    assert fit["success"] is True

    for key, true_val in ("a0", a0_true), ("a1", a1_true), ("b1", b1_true):
        rel = abs(fit[key] - true_val) / true_val
        assert rel <= 0.2

    Z_fit = Z_debye_k1(omega_pos, fit["a0"], fit["a1"], fit["b1"])
    R0_fit = (1.0 - Z_fit) / (1.0 + Z_fit)
    resid_R, resid_Z = passivity_residuals(omega_pos, Z_fit, R0_fit)
    assert resid_R <= 1e-12
    assert resid_Z <= 1e-12
