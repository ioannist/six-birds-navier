import cmath

from bhwave.barrier import barrier_reflection_transmission, poeschl_teller_V
from bhwave.echo import cavity_H, predicted_R_out_x0
from bhwave.scattering import scattering_reflection


def _setup_barrier():
    V0 = 1.0
    a = 1.0
    x_peak = 0.0

    def V(x: float) -> float:
        return poeschl_teller_V(x, V0, a, x_peak)

    return V


def _Z_debye(omega: float, a0: float, a1: float, b1: float) -> complex:
    return a0 + a1 / (1.0 - 1j * omega / b1)


def test_bh15_freqdep_reflection_formula():
    V = _setup_barrier()
    x0 = -8.0
    x_min = -12.0
    x_max = 12.0
    x_ref = 0.0
    Dx = x_ref - x0

    a0, a1, b1 = 0.2, 0.5, 2.0

    for omega in [0.5, 1.0, 2.0, 5.0]:
        Z = _Z_debye(omega, a0, a1, b1)
        R0 = (1.0 - Z) / (1.0 + Z)

        barrier = barrier_reflection_transmission(
            omega,
            V,
            x_min=x_min,
            x_max=x_max,
            x_ref=x_ref,
        )
        r = barrier.R_barrier
        t = barrier.T_barrier

        full = scattering_reflection(
            omega,
            V,
            x0,
            x_max,
            Z=Z,
        )
        R_full = full.R_out

        R_pred = predicted_R_out_x0(omega, R0, r, t, Dx)
        assert abs(R_full - R_pred) < 5e-6

        if abs(t) >= 0.2 and abs(R0) >= 0.05:
            H_pred = cavity_H(omega, R0, r, Dx)
            H_num = (R_full - r * cmath.exp(-2j * omega * Dx)) / (t * t * R0)
            rel_err = abs(H_num - H_pred) / max(1e-12, abs(H_pred))
            assert rel_err < 2e-5
