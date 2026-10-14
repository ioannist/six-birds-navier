#!/usr/bin/env python3
"""Exact Gaussian directional-derivative controls, not nonlinear evolution."""
from fractions import Fraction as F
from pathlib import Path
import json
from heat_ledger_controls import add, scale, multiply, derivative, dilate, gaussian_average


def poly_add(*polys):
    out = [F(0)] * max(map(len, polys))
    for p in polys:
        for j, c in enumerate(p):
            out[j] += c
    while len(out)>1 and not out[-1]:
        out.pop()
    return out


def poly_multiply(p, q):
    out = [F(0)] * (len(p)+len(q)-1)
    for i, c in enumerate(p):
        for j, d in enumerate(q):
            out[i+j] += c*d
    return out


def main():
    # Pi_beta/pi^(3/2)=n(beta)/d(beta)^(11/2).
    n=[F(0),F(192),F(640)]
    d=[F(3),F(16)]
    nprime=[j*c for j,c in enumerate(n)][1:]
    diff_numerator=poly_add(poly_multiply(nprime,d),[-F(11,2)*16*c for c in n])
    pure_heat_numerator=poly_add([-32*c for c in poly_multiply(n,d)],
                                 [F(0)]+[-4*c for c in diff_numerator])
    assert pure_heat_numerator == [F(0),F(-20736),F(-119808),F(-184320)]
    assert all(c<0 for c in pure_heat_numerator[1:])

    # A separate direct 3D Gaussian integration of the STATE derivative.
    p=[{(1,0,2):F(4),(1,0,0):F(-2)},
       {(0,1,2):F(-2),(0,1,0):F(1)},
       {(2,0,1):F(-4),(0,2,1):F(2),(0,0,1):F(1)}]
    b=[add(*(multiply(p[j],derivative(p[i],j)) for j in range(3))) for i in range(3)]
    lap=[add(*(derivative(derivative(pi,j),j) for j in range(3))) for pi in p]
    db=[add(*(add(multiply(lap[j],derivative(p[i],j)),
                   multiply(p[j],derivative(lap[i],j))) for j in range(3))) for i in range(3)]
    receiver_lambda=F(1,2)  # beta=3/8, H_(2beta)
    receiver=[scale(dilate(pi,receiver_lambda),receiver_lambda**6) for pi in p]
    receiver_lap=[add(*(derivative(derivative(pi,j,gaussian=receiver_lambda**2),
                                  j,gaussian=receiver_lambda**2) for j in range(3)))
                  for pi in receiver]
    q=2+receiver_lambda**2
    cubic=gaussian_average(add(*(add(multiply(receiver_lap[i],b[i]),
                                     multiply(receiver[i],db[i])) for i in range(3))),q)
    assert cubic == F(-53,729)
    # The terminal m2 edge is genuinely nonzero on the original datum.
    production=gaussian_average(add(*(multiply(lap[i],b[i]) for i in range(3))),3)
    assert 2*production == F(64,27)
    price_exponents={}
    for gamma in (F(3),F(7,2),F(4)):
        parabolic_price=(8*gamma-2*(4*gamma+1))/3
        assert parabolic_price == F(-2,3)
        price_exponents[str(gamma)]={
            'parabolic_cadence_native_price_exponent':str(parabolic_price),
            'zero_exponent_cadence_threshold':str(8*gamma/(4*gamma+1))}
    result={
        'scope':'Exact original Schwartz initial-state directional cubic controls; no nonlinear physical trajectory or terminal signed-source estimate.',
        'pure_heat_directional_flux_numerator_coefficients':list(map(str,pure_heat_numerator)),
        'pure_heat_directional_flux':'T_beta[Delta U]=-2304*pi^(3/2)*beta*(80*beta^2+52*beta+9)/(16*beta+3)^(13/2)',
        'pure_heat_directional_flux_sign':'Strictly negative at every beta>0; beta0 is zero.',
        'direct_beta_3_over_8_coefficient_of_(pi/(9/4))^(3/2)':str(cubic),
        'terminal_m2_scale_edge_coefficient_of_(pi/3)^(3/2)':str(-2*production),
        'true_physical_derivative':'At original initial datum bU: nu*b^3*T_beta[Delta U]+b^4*Q_beta(U), with full pressure quartic Q finite; small positive b gives nonzero negative derivative.',
        'uncomputed':'Full quartic pressure coefficient Q_beta and its terminal-uniform budget; no claim these are paid.',
        'lag_band_charge_algebraic_controls':price_exponents,
        'lag_band_charge_scope':'R^(-2/3) summability is a control of the derived estimate, not an actual Navier trajectory or disjoint-event construction.',
    }
    output=Path(__file__).resolve().parents[3]/'python/artifacts/unconditional-navier/heat-window-controls.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
