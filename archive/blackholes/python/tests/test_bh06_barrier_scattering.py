import numpy as np

from bhwave.barrier import barrier_reflection_transmission, poeschl_teller_V


def test_free_case_reflection_transmission():
    def V0(x: float) -> float:
        return 0.0

    x_min = -10.0
    x_max = 10.0
    x_ref = 0.0
    for omega in [0.7, 2.0, 4.5]:
        res = barrier_reflection_transmission(
            omega,
            V0,
            x_min=x_min,
            x_max=x_max,
            x_ref=x_ref,
        )
        assert res.success is True
        assert abs(res.R_barrier) < 1e-6
        assert abs(res.T_barrier - 1.0) < 1e-6
        assert abs(abs(res.R_barrier) ** 2 + abs(res.T_barrier) ** 2 - 1.0) < 1e-6


def test_poeschl_teller_unitarity_and_trend():
    V0 = 1.0
    a = 1.0
    x_peak = 0.0

    def V(x: float) -> float:
        return poeschl_teller_V(x, V0, a, x_peak)

    x_min = -12.0
    x_max = 12.0
    x_ref = 0.0
    omegas = [0.5, 1.0, 2.0, 5.0]
    results = {}
    for omega in omegas:
        res = barrier_reflection_transmission(
            omega,
            V,
            x_min=x_min,
            x_max=x_max,
            x_ref=x_ref,
        )
        assert res.success is True
        residual = abs(abs(res.R_barrier) ** 2 + abs(res.T_barrier) ** 2 - 1.0)
        assert residual < 1e-4
        assert abs(res.R_barrier) <= 1.0 + 1e-6
        results[omega] = res

    assert abs(results[5.0].R_barrier) < abs(results[0.5].R_barrier)
