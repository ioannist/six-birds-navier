"""Finite exact controls for analytic multipliers, Fourier units and symmetry."""
from pathlib import Path
import json
import sympy as s


def controls():
    r = s.symbols('r', nonnegative=True)
    derivative_factor = s.factor(s.diff(r*s.exp(-r*r+r), r)/s.exp(-r*r+r))
    assert s.simplify(derivative_factor+(r-1)*(2*r+1)) == 0
    assert (r*s.exp(-r*r+r)).subs(r, 1) == 1
    kappa, radius = s.symbols('kappa radius', positive=True)
    physical_ball = s.integrate(4*s.pi*r*r*(2*s.pi*r)**2, (r, 0, radius))
    angular_ball = s.simplify(physical_ball.subs(radius, kappa/(2*s.pi)))
    assert physical_ball == 16*s.pi**3*radius**5/5
    assert angular_ball == kappa**5/(10*s.pi**2)
    x, y, z = s.symbols('x y z', real=True)
    potential = x*s.exp(-x*x-y*y-z*z)
    velocity = s.Matrix([s.diff(potential, y), -s.diff(potential, x), 0])
    assert s.simplify(s.diff(velocity[0], x)+s.diff(velocity[1], y)) == 0
    origin = velocity.subs({x:0, y:0, z:0})
    assert origin == s.Matrix([0, -1, 0])
    half_turn = s.diag(-1, -1, 1)
    assert half_turn*origin != origin
    d, a, b, n = s.symbols('d a b n', nonnegative=True)
    quotient_bound = -2*d*d+2*d*(a*s.sqrt(n)+b)
    square_gap = s.expand((a*s.sqrt(n)+b)**2/2-quotient_bound)
    assert s.simplify(square_gap-2*(d-(a*s.sqrt(n)+b)/2)**2) == 0
    return {
        'status':'PASS',
        'scope':'Exact finite multiplier derivative, physical/angular Fourier volume conversion, solenoidal Schwartz initial jet and quotient Young identity. No NS evolution, functional fixed-point/exhaustion, backward uniqueness or Clay closure certified.',
        'heat_derivative_factor':str(derivative_factor),
        'physical_ball_coefficient':str(physical_ball),
        'angular_ball_coefficient':str(angular_ball),
        'nonsymmetric_origin':list(map(str,origin)),
        'quotient_square_gap':str(square_gap),
    }


if __name__ == '__main__':
    result=controls()
    artifact=Path(__file__).resolve().parents[2]/'artifacts/unconditional-navier/analytic-source-controls.json'
    artifact.parent.mkdir(parents=True,exist_ok=True)
    artifact.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
