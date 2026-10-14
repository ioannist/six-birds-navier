"""Exact initial jets for actual compact-data cofactor sign controls.

The mathematical compactification and local-time Taylor adapter are separate
proof duties, supplied in round58/parent-actual-cofactor-control.md.
"""
from pathlib import Path
import json
import sympy as s


def controls():
    x, y, z, a, b, nu = s.symbols("x y z a b nu", real=True, nonzero=True)
    coords = (x, y, z)
    position = s.Matrix(coords)
    origin = dict.fromkeys(coords, 0)

    def curl(v):
        return s.Matrix([s.diff(v[2], y)-s.diff(v[1], z),
                         s.diff(v[0], z)-s.diff(v[2], x),
                         s.diff(v[1], x)-s.diff(v[0], y)])

    def lap(v):
        return v.applyfunc(lambda q: sum(s.diff(q, c, 2) for c in coords))

    def zero(v):
        assert all(s.expand(q) == 0 for q in v)

    phi = -a*x**3/3+a*x*(y*y+z*z)/2
    quadratic = s.Matrix([s.diff(phi, c) for c in coords])
    strain = quadratic.jacobian(coords)
    strain_jets = [strain.diff(c) for c in coords]
    zero(curl(quadratic))
    assert s.expand(s.trace(strain)) == 0
    cofactor_sum = sum((q.adjugate() for q in strain_jets), s.zeros(3))
    assert cofactor_sum == s.diag(a*a, -3*a*a, -3*a*a)
    cases = []
    for axis, expected in ((0, 8), (1, -24)):
        omega0 = s.zeros(3, 1)
        omega0[axis] = b
        rotation = omega0.cross(position)/2
        matrix = s.Matrix.hstack(*(q*omega0 for q in strain_jets))
        zero(matrix-matrix.T)
        assert s.trace(matrix) == 0
        F = matrix*position
        radius2 = position.dot(position)
        quartic = -radius2*position.cross(F)/(42*nu)
        omega4 = (5*radius2*F-2*position*position.dot(F))/(42*nu)
        zero(curl(quartic)-omega4)
        zero(lap(omega4)-F/nu)
        u = quadratic+rotation+quartic
        omega = curl(u)
        zero(u.subs(origin))
        zero(omega.subs(origin)-omega0)
        zero(omega.jacobian(coords).subs(origin))
        jac = u.jacobian(coords)
        # Exact original vorticity equation at the initial face.
        omega_t = -omega.jacobian(coords)*u+jac*omega+nu*lap(omega)
        omega_gradient_t = omega_t.jacobian(coords).subs(origin)
        zero(omega_gradient_t-2*matrix)
        # Specific G starts at I: G_k,t=-2 S_k. Thus (G omega)_k,t=0.
        e_gradient_t = s.Matrix.hstack(*(-2*q*omega0 for q in strain_jets))+omega_gradient_t
        zero(e_gradient_t)
        for q in strain_jets:
            assert s.trace(q) == 0
        # The dissipative completed square is O(t^4); this is its t^2 source.
        leading = s.expand(8*nu*(omega0.T*cofactor_sum*omega0)[0])
        assert s.simplify(leading-expected*nu*a*a*b*b) == 0
        potentials = [-position.cross(v)/(degree+2)
                      for v, degree in ((quadratic, 2), (rotation, 1), (quartic, 4))]
        zero(curl(sum(potentials, s.zeros(3, 1)))-u)
        cases.append(dict(axis=axis+1, operator_r_over_detG_t2=str(leading),
                          omega_gradient_t=[str(q) for q in omega_gradient_t],
                          e_gradient_t_residual="0", divergence="0"))
    ta, tb, tc, te, tf, th, ti = s.symbols("ta tb tc te tf th ti", real=True)
    td, tg, tj = -ta-tf, -tb-ti, -tc-th
    E = ta**2+tg**2+tj**2+3*(tb**2+tc**2+td**2+tf**2+th**2+ti**2)+6*te**2
    N = tc**2+th**2+tj**2+2*(te**2+tf**2+ti**2)
    squares = 4*(ta+3*tf/4)**2+4*(tb+ti/4)**2+s.Rational(5,12)*(tf**2+ti**2)+s.Rational(2,3)*(tc-th)**2+s.Rational(8,3)*te**2
    assert s.expand(E-s.Rational(5,3)*N-squares) == 0
    harmonic_cases = []
    for phi_h, expected_E, expected_N in ((z**3/6-z*(x*x+y*y)/4, s.Rational(5,2), s.Rational(3,2)), (x*y*z, s.Integer(6), s.Integer(2))):
        jets = [s.hessian(phi_h, coords).diff(c) for c in coords]
        energy = sum(sum(entry**2 for entry in q) for q in jets)
        mass = sum(q[:, 2].dot(q[:, 2]) for q in jets)
        assert energy == expected_E and mass == expected_N
        harmonic_cases.append(dict(E=str(energy), N=str(mass), cofactor_t2=str(4*(mass-energy/2)), full_operator_t2_over_nu=str(2*(3*mass-2*energy))))
    return dict(status="PASS", scope="Exact original initial polynomial jets, vector-potential identities and restricted harmonic tensor cancellation; compactification and smooth local-time Taylor bridge are proved separately. No PDE simulation, singular solution or unconditional closure.",
                strain_cofactor_sum=[[str(q) for q in row] for row in cofactor_sum.tolist()], cases=cases,
                harmonic_sum_of_squares_residual="0", harmonic_cases=harmonic_cases)


if __name__ == "__main__":
    result = controls()
    target = Path(__file__).resolve().parents[2]/"artifacts/unconditional-navier/actual-cofactor-controls.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(dict(status=result["status"], scope=result["scope"], artifact=str(target)), indent=2))
