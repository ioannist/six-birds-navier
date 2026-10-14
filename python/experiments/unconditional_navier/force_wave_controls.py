"""Exact algebra controls; no singular PDE realization or proof import."""
from pathlib import Path
import json
import sympy as s


def controls():
    r, radius = s.symbols('r radius', positive=True)
    phi_inside = radius**2 / 2 - r**2 / 6
    phi_outside = radius**3 / (3*r)
    assert s.simplify(phi_inside.subs(r, radius)-phi_outside.subs(r, radius)) == 0
    assert s.simplify(s.diff(phi_inside, r).subs(r, radius)-s.diff(phi_outside, r).subs(r, radius)) == 0
    assert s.diff(phi_inside, r, 2) == -s.Rational(1, 3)
    assert s.simplify(s.diff(phi_outside, r, 2)-2*radius**3/(3*r**3)) == 0
    energy = s.Function('E')(r)
    integral_tail = s.Function('I')(r)
    # Given I'=-E/r^4, the head-derived tail differential is exact.
    derivative = s.diff(r**2*integral_tail, r).subs(s.diff(integral_tail, r), -energy/r**4)
    assert s.simplify(derivative-(2*r*integral_tail-energy/r**2)) == 0
    c0four = s.Rational(1024, 27)
    assert 4*s.Rational(27, 32)*c0four == 128
    # L2 volume changes: squared force amplitude B^-6 times volume B^3.
    assert -6+3 == -3
    aa = s.Matrix(3, 3, s.symbols('a0:9'))
    kk = s.Matrix(s.symbols('k0:3', real=True))
    amp = s.Matrix(s.symbols('v0:3', real=True))
    nu = s.symbols('nu', positive=True)
    k2 = kk.dot(kk)
    dotk = -aa.T*kk
    dota = -aa*amp+2*kk*(kk.dot(aa*amp))/k2-nu*k2*amp
    assert s.simplify(dotk.dot(amp)+kk.dot(dota)+nu*k2*kk.dot(amp)) == 0
    assert s.simplify(amp.dot(dota)+amp.dot(aa*amp)+nu*k2*amp.dot(amp)-2*amp.dot(kk)*kk.dot(aa*amp)/k2) == 0
    rr, zz = s.symbols('r z', positive=True)
    ur = s.Function('ur')(rr, zz)
    uz = s.Function('uz')(rr, zz)
    ut = s.Function('ut')(rr, zz)
    conservative = s.diff(rr**2*ur*ut, rr)/rr**2+s.diff(uz*ut, zz)
    advective = ur*s.diff(ut, rr)+uz*s.diff(ut, zz)+ur*ut/rr
    divergence = s.diff(ur, rr)+ur/rr+s.diff(uz, zz)
    assert s.simplify(conservative-advective-ut*divergence) == 0
    return {
        'status': 'PASS',
        'scope': 'Exact Newton-ball potential, radial-tail derivative, coarse forced Young coefficient, transverse affine-wave identity, and cylindrical divergence algebra. No primary theorem proof, actual forced-profile asymptotics, source supplier, or Clay closure.',
        'forced_record_Q_coefficient': 128,
        'force_global_L2_exponent': '-3/2',
        'newton_ball_inside_head_coefficient': '1/6',
        'affine_wave': 'transversality and amplitude work, with k nonzero',
    }


if __name__ == '__main__':
    result = controls()
    artifact = Path(__file__).resolve().parents[2]/'artifacts/unconditional-navier/force-wave-controls.json'
    artifact.parent.mkdir(parents=True, exist_ok=True)
    artifact.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
