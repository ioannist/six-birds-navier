"""Exact Gaussian and temporal endpoint controls, without a nonlinear PDE claim."""
from pathlib import Path
import json
import sympy as s


def controls():
    x, y, z = s.symbols('x y z', real=True)
    radius2 = x*x+y*y+z*z
    h = s.Matrix([(1-2*y*y)*s.exp(-radius2),
                  2*x*y*s.exp(-radius2), 0])
    assert s.simplify(s.diff(h[0], x)+s.diff(h[1], y)) == 0
    # Gaussian polynomial expectation, covariance (2a)^-1 I.
    a = s.symbols('a', positive=True)
    polynomial = s.expand(h.dot(h)*s.exp(2*radius2))
    def gaussian_integral(poly, rate):
        result = 0
        for powers, coefficient in s.Poly(poly, x, y, z).terms():
            product = coefficient
            for power in powers:
                if power % 2:
                    product = 0
                    break
                product *= s.gamma(s.Rational(power+1, 2))/rate**s.Rational(power+1, 2)
            result += product
        return s.simplify(result)
    norm2 = gaussian_integral(polynomial, s.Integer(2))
    assert s.simplify(norm2-(s.pi/2)**s.Rational(3, 2)) == 0
    heat0 = gaussian_integral(1-2*y*y, s.Rational(5, 4))/(4*s.pi)**s.Rational(3, 2)
    assert s.simplify(heat0-5**(-s.Rational(5, 2))) == 0
    assert gaussian_integral(2*x*y, s.Rational(5, 4)) == 0
    ell = s.symbols('ell', positive=True)
    temporal_mass = s.integrate(ell**-2, (ell, s.Rational(3, 2), s.oo))
    assert temporal_mass == s.Rational(2, 3)
    n = s.symbols('n', positive=True)
    endpoint = s.integrate(ell**(-s.Rational(1, 2)), (ell, 2, n))
    assert s.simplify(endpoint-2*(s.sqrt(n)-s.sqrt(2))) == 0
    j = s.symbols('j', positive=True)
    shell_mass = s.integrate(ell**-2, (ell, j, j+1))
    assert s.simplify(shell_mass-1/(j*(j+1))) == 0
    return {
        'status': 'PASS',
        'scope': 'Exact solenoidal Gaussian norm/heat value and logarithmic temporal integrals. No actual Reynolds-force realization, nonlinear NS evolution, endpoint source supplier or Clay closure certified.',
        'gaussian_norm_squared': str(norm2),
        'unit_heat_origin_first_component': str(heat0),
        'uniform_fourth_power_temporal_bound': str(temporal_mass),
        'terminal_positive_response_factor': str(endpoint),
        'logarithmic_shell_mass': str(shell_mass),
    }


if __name__ == '__main__':
    result = controls()
    artifact = Path(__file__).resolve().parents[2]/'artifacts/unconditional-navier/mean-force-endpoint-controls.json'
    artifact.parent.mkdir(parents=True, exist_ok=True)
    artifact.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
