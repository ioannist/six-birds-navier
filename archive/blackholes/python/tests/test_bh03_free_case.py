from bhwave.conventions import R_from_Z
from bhwave.scattering import scattering_reflection


def V0(x: float) -> float:
    return 0.0


def test_free_case_reflection_matches_boundary():
    x0 = 0.0
    x_max_list = [3.0, 7.0]
    omega_list = [0.7, 2.0, 4.5]
    Z_list = [1.0 + 0j, 0.0 + 0j, 0.2 + 0.3j, 0.0 + 0.5j]

    for omega in omega_list:
        for x_max in x_max_list:
            for Z in Z_list:
                res = scattering_reflection(
                    omega,
                    V0,
                    x0,
                    x_max,
                    Z=Z,
                )
                assert res.success is True
                R_expected = R_from_Z(Z)
                assert abs(res.R_out - R_expected) < 1e-6
                if Z.real == 0.0 and Z.imag != 0.0:
                    assert abs(abs(res.R_out) - 1.0) < 1e-6


def test_free_case_independence_from_xmax():
    x0 = 0.0
    omega = 2.0
    Z = 0.2 + 0.3j

    res_a = scattering_reflection(omega, V0, x0, 3.0, Z=Z)
    res_b = scattering_reflection(omega, V0, x0, 7.0, Z=Z)

    assert res_a.success is True
    assert res_b.success is True
    assert abs(res_a.R_out - res_b.R_out) < 1e-7
