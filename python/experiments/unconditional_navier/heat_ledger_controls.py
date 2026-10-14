#!/usr/bin/env python3
"""Exact polynomial/Gaussian initial-state controls; no simulated NS path."""
from fractions import Fraction as F
from math import comb, isqrt
from pathlib import Path
import json

ZERO = (0, 0, 0)


def add(*polys):
    out = {}
    for poly in polys:
        for n, c in poly.items():
            out[n] = out.get(n, F(0)) + c
    return {n: c for n, c in out.items() if c}


def scale(poly, c):
    return {n: c * v for n, v in poly.items() if c * v}


def multiply(a, b):
    out = {}
    for n, c in a.items():
        for m, d in b.items():
            k = tuple(x + y for x, y in zip(n, m))
            out[k] = out.get(k, F(0)) + c * d
    return {n: c for n, c in out.items() if c}


def derivative(poly, j, gaussian=1):
    terms = []
    for n, c in poly.items():
        if n[j]:
            k = list(n); k[j] -= 1
            terms.append({tuple(k): c * n[j]})
        k = list(n); k[j] += 1
        terms.append({tuple(k): -2 * gaussian * c})
    return add(*terms)


def double_factorial(n):
    out = 1
    while n > 0:
        out *= n; n -= 2
    return out


def gaussian_average(poly, q):
    """Integral divided by (pi/q)^(3/2), with exact positive rational q."""
    return sum((c * F(double_factorial(n[0]-1) * double_factorial(n[1]-1)
                         * double_factorial(n[2]-1), (2*q)**(sum(n)//2))
                for n, c in poly.items() if all(k % 2 == 0 for k in n)), F(0))


def linear_power(constant, degree):
    return [F(comb(degree, j)) * constant**(degree-j) for j in range(degree+1)]


def polynomial_value(poly, t):
    return sum((c*t**j for j, c in enumerate(poly)), F(0))


def dilate(poly, factor):
    return {n: c*factor**sum(n) for n, c in poly.items()}


def factor_linear(poly, root):
    """Exact division by t-root; only called after checking that root."""
    assert polynomial_value(poly, root) == 0
    result = [F(0)] * (len(poly)-1)
    result[-1] = poly[-1]
    for j in range(len(poly)-2, 0, -1):
        result[j-1] = poly[j] + root*result[j]
    assert poly[0] + root*result[0] == 0
    return result


def main():
    # U_i=p_i exp(-r^2), exactly the historical real Schwartz control.
    p = [{(1, 0, 2): F(4), (1, 0, 0): F(-2)},
         {(0, 1, 2): F(-2), (0, 1, 0): F(1)},
         {(2, 0, 1): F(-4), (0, 2, 1): F(2), (0, 0, 1): F(1)}]
    assert not add(*(derivative(p[j], j) for j in range(3)))
    b = [add(*(multiply(p[j], derivative(p[i], j)) for j in range(3)))
         for i in range(3)]
    lap = [add(*(derivative(derivative(pi, j), j) for j in range(3))) for pi in p]
    for i in range(3):
        dilation = add(scale(p[i], -12), *(scale(multiply(
            {tuple(int(k == j) for k in range(3)): F(1)},
            derivative(p[i], j)), -2) for j in range(3)))
        assert lap[i] == dilation  # heat generator of lambda^6 U(lambda x)
    assert gaussian_average(add(*(multiply(p[i], b[i]) for i in range(3))), 3) == 0
    production = gaussian_average(add(*(multiply(lap[i], b[i]) for i in range(3))), 3)
    assert production == F(32, 27)

    # t=lambda^2=(1+8theta)^(-1); receiver has lambda^7 or lambda^9.
    terms = []
    for i in range(3):
        for n, c in p[i].items():
            receiver_power = (sum(n)-1)//2
            assert sum(n) in (1, 3)
            for k, d in b[i].items():
                degree = tuple(a+b for a, b in zip(n, k))
                if any(v % 2 for v in degree):
                    continue
                moment = sum(degree)//2
                coeff = c*d*F(double_factorial(degree[0]-1)
                              * double_factorial(degree[1]-1)
                              * double_factorial(degree[2]-1), 2**moment)
                terms.append((receiver_power, moment, coeff))
    denominator_power = max(moment for _, moment, _ in terms)
    numerator = [F(0)] * (denominator_power+2)
    for shift, moment, coeff in terms:
        for j, factor in enumerate(linear_power(2, denominator_power-moment)):
            numerator[j+shift] += coeff*factor
    while len(numerator) > 1 and not numerator[-1]:
        numerator.pop()
    assert polynomial_value(numerator, F(1)) == 0
    initial_derivative = -8*sum((j*c for j, c in enumerate(numerator)), F(0))/3**denominator_power
    assert initial_derivative == 2*production
    while polynomial_value(numerator, F(-2)) == 0:
        numerator = factor_linear(numerator, F(-2)); denominator_power -= 1
    # N(t)=(1-t) times a strictly positive Bernstein polynomial on [0,1].
    positive = [-c for c in factor_linear(numerator, F(1))]
    d = len(positive)-1
    bernstein = [sum((positive[j]*F(comb(k, j), comb(d, j))
                     for j in range(k+1)), F(0)) for k in range(d+1)]
    assert all(c > 0 for c in bernstein)
    energy = gaussian_average(add(*(multiply(pi, pi) for pi in p)), 2)
    z0 = gaussian_average(add(*(multiply(derivative(pi, j), derivative(pi, j))
                               for pi in p for j in range(3))), 2)
    assert z0 == F(117, 8)
    # Nonzero actual H^3 nonlinear work; this tests the full alternating
    # heat-source readout via its separately proved epsilon-to-zero bridge.
    third_power = p
    for _ in range(3):
        third_power = [scale(add(*(derivative(derivative(pi, j), j)
                                   for j in range(3))), -1) for pi in third_power]
    h3_work = gaussian_average(add(*(multiply(third_power[i], b[i])
                                    for i in range(3))), 3)
    assert h3_work == F(-19712, 27)
    # Direct Gaussian integration of actual rescaled fields, not just an
    # exponent string. Square R keeps the fixed-energy amplitude rational.
    receiver_lambda = F(1, 2)  # c=3/8 for H_(2c)
    qbase = 2 + receiver_lambda**2
    receiver = [scale(dilate(pi, receiver_lambda), receiver_lambda**6) for pi in p]
    base_flux = gaussian_average(add(*(multiply(receiver[i], b[i])
                                      for i in range(3))), qbase)
    assert base_flux > 0
    for radius in (4, 9):
        amplitude = radius*isqrt(radius)
        scaled_field = [scale(dilate(pi, radius), amplitude) for pi in p]
        scaled_b = [add(*(multiply(scaled_field[j], derivative(
            scaled_field[i], j, gaussian=radius**2)) for j in range(3)))
                    for i in range(3)]
        scaled_receiver = [scale(dilate(pi, radius*receiver_lambda),
                                  amplitude*receiver_lambda**6) for pi in p]
        # Normalize all spatial integrals by the same base Gaussian factor.
        scaled_flux = gaussian_average(add(*(multiply(scaled_receiver[i], scaled_b[i])
                                            for i in range(3))), radius**2*qbase)/radius**3
        scaled_energy = gaussian_average(add(*(multiply(pi, pi) for pi in scaled_field)),
                                         2*radius**2)/radius**3
        scaled_z = gaussian_average(add(*(multiply(
            derivative(pi, j, gaussian=radius**2),
            derivative(pi, j, gaussian=radius**2))
            for pi in scaled_field for j in range(3))), 2*radius**2)/radius**3
        assert scaled_flux == radius**2*isqrt(radius)*base_flux
        assert scaled_energy == energy
        assert scaled_z == radius**2*z0
    result = {
        "scope": "Exact original Schwartz initial-state control; no numerical evolution, fixed-data terminal singularity or all-data source budget.",
        "heat_response": "H_theta U=lambda^6 U(lambda*x), lambda=(1+4theta)^(-1/2)",
        "t_parameter": "t=(1+8theta)^(-1)",
        "flux_formula": "Pi_theta=pi^(3/2)*t^(7/2)*N(t)/(2+t)^(3/2+d)",
        "d": denominator_power,
        "N_coefficients_ascending": list(map(str, numerator)),
        "positive_factor_coefficients_ascending": list(map(str, positive)),
        "positive_factor_Bernstein_coefficients": list(map(str, bernstein)),
        "strict_flux_sign": "Pi_theta(U)>0 for every theta>0; Pi_theta(-U)=-Pi_theta(U)",
        "E0_coefficient_of_(pi/2)^(3/2)": str(energy),
        "Z0_coefficient_of_(pi/2)^(3/2)": str(z0),
        "Pi_derivative_at_zero_coefficient_of_(pi/3)^(3/2)": str(2*production),
        "pairing_(-Delta)^3_U_with_B_coefficient_of_(pi/3)^(3/2)": str(h3_work),
        "fixed_energy_scale": "U_R=R^(3/2)U(R*x); theta_R=c/R^2; Pi=R^(5/2)Pi_c, Z-Ztheta=R^2 Z0*(1-(1+4c)^(-11/2)), E-Etheta=E0*(1-(1+4c)^(-9/2))",
        "direct_scaled_Gaussian_controls": "Exact spatial integrations at R=4,9 and c=3/8 verify Pi scaling, fixed energy and enstrophy scaling; the all-R bridge is written in source-gate.md.",
    }
    root = Path(__file__).resolve().parents[3]
    output = root / 'python/artifacts/unconditional-navier/heat-ledger-controls.json'
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
