"""Exact source-budget and Bernoulli-weight controls, not an NS simulation."""
from fractions import Fraction as F
from pathlib import Path
import json
import sympy as s


def controls():
    r = s.symbols("r", positive=True)
    # Actual window Young maximum; C0^4=1024/27 belongs to the cubic law.
    c0four = s.Rational(1024, 27)
    assert s.Rational(27, 16)*c0four == 64
    theta = 1/s.sqrt(1+r*r)
    A = (r*s.sqrt(1+r*r)-s.asinh(r))/2
    assert s.simplify(s.diff(A, r)-r*r*theta) == 0
    # These exact monotonicity derivatives prove the integrated comparisons.
    assert s.simplify(s.diff(theta, r)+r/(1+r*r)**s.Rational(3, 2)) == 0
    assert s.simplify(s.diff(r*theta, r)-(1+r*r)**(-s.Rational(3, 2))) == 0
    at_zero = s.limit(theta/2-A/r**3, r, 0)
    ar_zero = s.limit(2*A/r**3-theta/2, r, 0)
    assert at_zero == ar_zero == s.Rational(1, 6)
    powers = []
    for alpha in (F(1, 2), F(1), F(2), F(5, 2)):
        powers.append(dict(alpha=str(alpha), tangential_over_theta=str((1-alpha)/(2*(3-alpha))),
                           radial_over_theta=str((1+alpha)/(2*(3-alpha))),
                           native_L6_weight_pairing_integrable=alpha>2))
    # Summable translated bump prices with nonsummable L^(9/2) price.
    amp, radius, center = F(-3), F(5), F(12)
    assert 2*amp+radius == -1
    assert 6*amp+3*radius == -3
    assert 2*amp+3*radius-center == -3
    assert F(9, 2)*amp+3*radius == F(3, 2)
    # Annular reverse norm control: beta=3/4 gives L9/2 but unpaid weight.
    beta = F(3, 4)
    assert 1-2*beta == F(-1, 2)
    assert 3-6*beta == F(-3, 2)
    assert 3-F(9, 2)*beta == F(-3, 8)
    assert 2-2*beta == F(1, 2)
    return dict(status="PASS", scope="Exact Gaussian-independent source coefficient, radial-weight calculus and disjoint-profile norm powers; no stationary NS realization, endpoint source, measure construction or unconditional closure.",
                arbitrary_shifted_window_Q_bound="64*length/nu^4+1/nu",
                radial_A=str(A), origin_head_eigenvalues=[str(at_zero), str(ar_zero)],
                formal_power_weights=powers,
                translated_profile_base2_exponents=dict(Z=-1, L6_sixth=-3, weighted_energy=-3, L9_2_power="3/2"),
                annular_profile_radius_exponents=dict(Z="-1/2", L6_sixth="-3/2", L9_2_power="-3/8", weighted_energy="1/2"))


if __name__ == "__main__":
    result = controls()
    target = Path(__file__).resolve().parents[2]/"artifacts/unconditional-navier/head-source-controls.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(dict(status=result["status"], scope=result["scope"], artifact=str(target)), indent=2))
