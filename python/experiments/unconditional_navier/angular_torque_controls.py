"""Exact radial cutoff and original initial-torque controls, not time evolution."""
from pathlib import Path
import json
import sympy as s


def controls():
    r, length = s.symbols('r length', positive=True)
    weight = r*r/(1+r*r/length**2)
    laplacian = s.diff(weight, r, 2)+2*s.diff(weight, r)/r
    expected = (6-2*r*r/length**2)/(1+r*r/length**2)**3
    assert s.simplify(laplacian-expected) == 0
    assert s.simplify(s.diff(weight, r)**2/(4*weight)-(1+r*r/length**2)**-3) == 0
    # The squared gradient ratio is <=1; Delta w <=6 follows algebraically.
    x = s.symbols('x', nonnegative=True)
    assert s.expand(6*(1+x)**3-(6-2*x)) == 6*x**3+18*x**2+20*x
    z, eps, kappa, theta = s.symbols('z eps kappa theta', real=True)
    gaussian = s.exp(-r*r-z*z)
    psi = eps*r*gaussian*(s.cos(theta)+kappa*r*r*s.sin(theta))
    ur, ut = s.diff(psi, theta)/r, -s.diff(psi, r)
    covariance = s.simplify(s.integrate(ur*ut, (theta, 0, 2*s.pi))/(2*s.pi))
    torque = s.simplify(-s.diff(r*r*covariance, r)/r)
    assert s.simplify(covariance-eps**2*kappa*r*r*gaussian**2) == 0
    assert s.simplify(torque+4*eps**2*kappa*r*r*(1-r*r)*gaussian**2) == 0
    assert s.integrate(r**3*(1-r*r)*s.exp(-2*r*r), (r, 0, s.oo)) == 0
    return {
        'status': 'PASS',
        'scope': 'Exact bounded radial weight derivative/gradient/Laplacian and mean-zero Schwartz initial covariance/torque algebra. No evolved moment persistence, endpoint source estimate, globality, or numerical NS trajectory certified.',
        'cutoff_laplacian': str(expected),
        'gradient_squared_ratio': '(1+r^2/L^2)^(-3) <= 1',
        'initial_reynolds_torque': str(torque),
        'signed_initial_torque_integral': 0,
    }


if __name__ == '__main__':
    result = controls()
    artifact = Path(__file__).resolve().parents[2]/'artifacts/unconditional-navier/angular-torque-controls.json'
    artifact.parent.mkdir(parents=True, exist_ok=True)
    artifact.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
