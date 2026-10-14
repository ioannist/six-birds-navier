#!/usr/bin/env python3
"""Exact original initial-jet controls; no simulated solution or source bound."""
from fractions import Fraction as F
from pathlib import Path
import json
import re
from heat_ledger_controls import add, scale, multiply, derivative

ORIGIN = (0, 0, 0)
POTENTIAL = (
    'x^4+x^3y+x^3z+x^3-x^2y+x^2z-x^2-xy^2z-xy^2+xyz^2-xz^2+xz-x'
    '+y^2z^2-y^2+yz^3-yz^2-yz+z^4-z^3-z^2+z+1',
    '-x^3+x^2y-x^2z^2+xy^2z-xy^2-xyz^2-xyz-xy-xz^3+y^4-y^3z-y^3'
    '-y^2z^2+y^2z-yz^3-yz-y+z^4-z^2-1',
    'x^4-x^3y-x^3z+x^3-x^2y-x^2z^2-x^2-xy^3-xyz^2+xz^3-xz^2-xz+x'
    '+y^4-y^3+y^2z^2-y^2z-y^2+yz^3-yz^2-yz+z^3-z',
)


def parse_unit_polynomial(text):
    """Parse only signed monomials of coefficient one, without evaluation."""
    pieces = []
    for term in re.findall(r'[+-]?[^+-]+', text):
        sign = -1 if term.startswith('-') else 1
        monomial = term.lstrip('+-')
        powers = [0, 0, 0]
        if monomial != '1':
            tokens = list(re.finditer(r'([xyz])(?:\^(\d+))?', monomial))
            assert ''.join(m.group(0) for m in tokens) == monomial
            for token in tokens:
                powers['xyz'.index(token.group(1))] += int(token.group(2) or 1)
        pieces.append({tuple(powers): F(sign)})
    return add(*pieces)


def d(p, j):
    return derivative(p, j, gaussian=0)


def curl(v, gaussian=0):
    return [add(derivative(v[(i+2) % 3], (i+1) % 3, gaussian),
                scale(derivative(v[(i+1) % 3], (i+2) % 3, gaussian), -1))
            for i in range(3)]


def lap(p):
    return add(*(d(d(p, j), j) for j in range(3)))


def dot(a, b):
    return sum((x*y for x, y in zip(a, b)), F(0))


def cross(a, b):
    return [a[(i+1) % 3]*b[(i+2) % 3]-a[(i+2) % 3]*b[(i+1) % 3]
            for i in range(3)]


def poly_cross(a, b):
    return [add(multiply(a[(i+1) % 3], b[(i+2) % 3]),
                scale(multiply(a[(i+2) % 3], b[(i+1) % 3]), -1))
            for i in range(3)]


def val(p, point=ORIGIN):
    return sum((c*point[0]**n[0]*point[1]**n[1]*point[2]**n[2]
                for n, c in p.items()), F(0))


def mv(a, v):
    return [dot(row, v) for row in a]


def controls():
    p = [parse_unit_polynomial(s) for s in POTENTIAL]
    u = curl(p)
    assert not add(*(d(u[i], i) for i in range(3)))
    w = curl(u)
    ap = [[d(u[i], j) for j in range(3)] for i in range(3)]
    vp = [[d(w[i], k) for i in range(3)] for k in range(3)]
    jp = [poly_cross(w, v) for v in vp]
    kp = add(*(multiply(t, t) for v in jp for t in v))
    a = [[val(t) for t in row] for row in ap]
    ww = [val(t) for t in w]
    v = [[val(t) for t in row] for row in vp]
    jj = [[val(t) for t in row] for row in jp]
    ta = mix = ga = vc = gj = F(0)
    for k in range(3):
        deform = [x+y for x, y in zip(cross(mv(a, ww), v[k]),
                                      cross(ww, mv(a, v[k])))]
        ta += dot(jj[k], deform)
        reduced = -dot(jj[k], mv(list(zip(*a)), jj[k]))
        assert dot(jj[k], deform) == reduced
        dak = [[val(d(t, k)) for t in row] for row in ap]
        ga += dot(jj[k], cross(ww, mv(dak, ww)))
        for ell in range(3):
            mix -= a[ell][k]*dot(jj[k], jj[ell])
        for j in range(3):
            djv = [val(d(t, j)) for t in vp[k]]
            vc += dot(jj[k], cross(v[j], djv))
            dj = [val(d(t, j)) for t in jp[k]]
            gj += dot(dj, dj)

    # Independently prolong D w=A w+nu Delta w before forming D K.
    # This checks the complete balance against direct PDE differentiation.
    awp = [add(*(multiply(ap[i][j], w[j]) for j in range(3)))
           for i in range(3)]
    aw = [val(t) for t in awp]
    lw = [val(lap(t)) for t in w]
    direct_inviscid = direct_viscous = F(0)
    for k in range(3):
        dv0 = [val(d(awp[i], k))-sum(a[ell][k]*v[ell][i]
                                   for ell in range(3)) for i in range(3)]
        dv1 = [val(lap(t)) for t in vp[k]]
        dj0 = [x+y for x, y in zip(cross(aw, v[k]), cross(ww, dv0))]
        dj1 = [x+y for x, y in zip(cross(lw, v[k]), cross(ww, dv1))]
        direct_inviscid += 2*dot(jj[k], dj0)
        direct_viscous += 2*dot(jj[k], dj1)
    delta_k = val(lap(kp))
    lhs = [direct_inviscid, direct_viscous-delta_k]
    rhs = [2*(ta+mix+ga), -4*vc-2*gj]
    assert lhs == rhs == [F(18940), F(-145124)]
    assert ga == 812 and vc == 15974 and val(kp) == 3558
    assert ww == [2, 1, 4] and v[2] == [3, -2, 0]

    # Actual Gaussian U at (1,2,0): every first derivative carries e^-5.
    ug = [{(1, 0, 2): F(4), (1, 0, 0): F(-2)},
          {(0, 1, 2): F(-2), (0, 1, 0): F(1)},
          {(2, 0, 1): F(-4), (0, 2, 1): F(2), (0, 0, 1): F(1)}]
    wg = [val(t, (1, 2, 0)) for t in curl(ug, gaussian=1)]
    ag = [[val(derivative(ug[i], j), (1, 2, 0)) for j in range(3)]
          for i in range(3)]
    omega2, stretch = dot(wg, wg), dot(wg, mv(ag, wg))
    assert omega2 == 144 and stretch == 720 and stretch/omega2 == 5
    threshold = F(1)/(16*F(16, 3)**3)
    assert threshold == F(27, 65536)
    # A second, simple solenoidal shear jet falsifies pointwise negativity
    # of the COMPLETE viscous bending source, rather than one isolated term.
    shear = [{}, {}, {(0, 1, 0): F(1), (2, 0, 0): F(-1, 2),
                     (2, 1, 0): F(1, 4)}]
    sw = curl(shear)
    sv = [[d(t, k) for t in sw] for k in range(3)]
    sj = [poly_cross(sw, t) for t in sv]
    sk = add(*(multiply(t, t) for row in sj for t in row))
    assert not add(*(d(shear[i], i) for i in range(3)))
    saw = [add(*(multiply(d(shear[i], j), sw[j]) for j in range(3)))
           for i in range(3)]
    assert not any(saw)  # Exact (omega.grad)u=0 on this local shear.
    sww, slw = [val(t) for t in sw], [val(lap(t)) for t in sw]
    scross = sgrad = sdirect = F(0)
    for k in range(3):
        sjk, svk = [val(t) for t in sj[k]], [val(t) for t in sv[k]]
        direct = [x+y for x, y in zip(cross(slw, svk),
                                     cross(sww, [val(lap(t)) for t in sv[k]]))]
        sdirect += 2*dot(sjk, direct)
        for j in range(3):
            dj = [val(d(t, j)) for t in sj[k]]
            sgrad += dot(dj, dj)
            scross += dot(sjk, cross([val(t) for t in sv[j]],
                                    [val(d(t, j)) for t in sv[k]]))
    assert val(sk) == 1 and val(lap(sk)) == 0
    assert sgrad == F(1, 2) and scross == F(-1, 2)
    assert sdirect-val(lap(sk)) == -2*sgrad-4*scross == 1
    return {
        'scope': 'Exact original initial jets after cutoff of potential; no time-integrated source or simulated path.',
        'polynomial_potential': list(POTENTIAL),
        'omega_origin': [str(t) for t in ww],
        'domega_columns': [[str(t) for t in row] for row in v],
        'K_origin': str(val(kp)),
        'contractions': dict(zip(('T_A', 'T_mix', 'T_gradA', 'V_cross', 'gradJ_squared'),
                                map(str, (ta, mix, ga, vc, gj)))),
        'complete_(D-nu_Delta)K_coefficients_constant_nu': list(map(str, lhs)),
        'direct_DK_coefficients_constant_nu': list(map(str, (direct_inviscid, direct_viscous))),
        'DeltaK_origin': str(delta_k),
        'Gaussian_point': '(1,2,0)',
        'Gaussian_omega2': '144 exp(-10)',
        'Gaussian_stretching': '720 exp(-15)',
        'Gaussian_alpha': '5 exp(-5), reverses for -U',
        'small_product_threshold_coefficient_of_nu4': str(threshold),
        'positive_viscous_bending_jet': {
            'local_velocity': '(0,0,y-x^2/2+x^2*y/4)',
            'omega': '(1+x^2/4,x-x*y/2,0)',
            'K_origin': '1', 'DeltaK_origin': '0',
            'gradJ_squared_origin': str(sgrad), 'V_cross_origin': str(scross),
            'complete_viscous_(D-nu_Delta)K_origin': '+nu',
            'scope': 'Falsifies pointwise nonpositivity only; cutoff of a curl potential preserves the actual Schwartz jet. No global integral or singularity assertion.',
        },
        'correctness_bridge': 'Plain polynomial derivatives, full direct material prolongation and independent split RHS agree; cutoff potential preserves these original solenoidal jets.',
        'unpaid': 'General local-direction dose, signed bending propagation, swirl source/shear doses, and unconditional original Clay(A).',
    }


if __name__ == '__main__':
    result = controls()
    output = Path(__file__).resolve().parents[2]/'artifacts/unconditional-navier/geometry-controls.json'
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
