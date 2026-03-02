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


def test_reflection_formula_match():
    V = _setup_barrier()
    x0 = -8.0
    x_min = -12.0
    x_max = 12.0
    x_ref = 0.0
    Dx = x_ref - x0
    R0 = 0.2 * cmath.exp(0.3j)

    for omega in [0.5, 1.0, 2.0, 5.0]:
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
            R_boundary=R0,
        )
        R_full = full.R_out

        R_pred = predicted_R_out_x0(omega, R0, r, t, Dx)
        assert abs(R_full - R_pred) < 2e-6


def test_H_factor_match():
    V = _setup_barrier()
    x0 = -8.0
    x_min = -12.0
    x_max = 12.0
    x_ref = 0.0
    Dx = x_ref - x0
    R0 = 0.2 * cmath.exp(0.3j)

    for omega in [1.0, 2.0]:
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
            R_boundary=R0,
        )
        R_full = full.R_out

        H_pred = cavity_H(omega, R0, r, Dx)
        H_num = (R_full - r * cmath.exp(-2j * omega * Dx)) / (t * t * R0)
        rel_err = abs(H_num - H_pred) / max(1e-12, abs(H_pred))
        assert rel_err < 1e-5
