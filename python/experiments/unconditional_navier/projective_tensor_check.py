"""Exact algebra controls for the parent projective tensor dossier.

These checks verify finite-dimensional Hessian/curl identities, not an
endpoint PDE estimate or the unconditional Clay claim.
"""
import sympy as sp


def main():
    x, y, z = sp.symbols("x y z", real=True)
    ell, kappa = sp.symbols("ell kappa", positive=True)
    radial, p, q = sp.symbols("radial p q", real=True)
    v = sp.Matrix([x, y, z])
    h = sp.Matrix([radial, p, q])
    tangent = sp.Matrix([0, p, q])
    axis = sp.Matrix([1, 0, 0])
    projection = axis * axis.T
    tensor = v * v.T / (sp.sqrt(v.dot(v)) + kappa)
    expected = (
        2*kappa**2/(ell+kappa)**3 * projection * radial**2
        + 2*kappa/(ell+kappa)**2
        * (axis*tangent.T + tangent*axis.T) * radial
        + 2/(ell+kappa) * tangent*tangent.T
        - ell/(ell+kappa)**2 * projection * tangent.dot(tangent)
    )
    residues = []
    for i in range(3):
        for j in range(3):
            actual = (h.T * sp.hessian(tensor[i, j], (x, y, z)) * h)[0]
            actual = actual.subs({x: ell, y: 0, z: 0})
            residues.append(sp.simplify(actual - expected[i, j]))
    assert all(value == 0 for value in residues), residues
    f, g = sp.Function("f")(x), sp.Function("g")(x)
    potential = sp.Matrix([f*y, 0, g])

    def curl(field):
        return sp.Matrix([
            sp.diff(field[2], y)-sp.diff(field[1], z),
            sp.diff(field[0], z)-sp.diff(field[2], x),
            sp.diff(field[1], x)-sp.diff(field[0], y),
        ])

    velocity = curl(potential)
    vorticity = curl(velocity)
    assert velocity == sp.Matrix([0, -sp.diff(g, x), -f])
    assert vorticity == sp.Matrix([0, sp.diff(f, x), -sp.diff(g, x, 2)])
    assert sum(sp.diff(velocity[i], c) for i, c in enumerate((x, y, z))) == 0
    print("PASS: nine exact contracted-Hessian residues; local curl/curl and div-free controls")
    print("Scope: algebra only; no endpoint source or Clay proof")


if __name__ == "__main__":
    main()
