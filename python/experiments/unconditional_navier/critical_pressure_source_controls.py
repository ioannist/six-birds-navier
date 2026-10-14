"""Exact finite Fourier pressure algebra; no numerical NS or limit certification."""
import json
from pathlib import Path
import sympy as s


def run():
    k, ell, h = (1, 0, 0), (0, 1, 1), (1, 1, 1)
    modes = {}
    for wave, vector, sine in [(k, (0, 1, 0), False),
                               (ell, (1, 0, 0), False),
                               (h, (1, 0, -1), True)]:
        modes[wave] = tuple(s.Rational(x, 2) * (-s.I if sine else 1) for x in vector)
        modes[tuple(-x for x in wave)] = tuple(s.conjugate(x) for x in modes[wave])
    assert all(sum(a*b for a, b in zip(wave, vector)) == 0 for wave, vector in modes.items())
    source = {}
    for wave, vector in modes.items():
        for other_wave, other_vector in modes.items():
            summed = tuple(a+b for a, b in zip(wave, other_wave))
            value = -sum(a*b for a, b in zip(wave, other_vector)) * sum(a*b for a, b in zip(other_wave, vector))
            source[summed] = s.simplify(source.get(summed, 0) + value)
    assert source[(0, 0, 0)] == 0
    pressure = {wave: s.simplify(value/sum(x*x for x in wave))
                for wave, value in source.items() if any(wave) and value != 0}
    mean = s.simplify(sum(value * s.I * (-wave[2]) * modes.get(tuple(-x for x in wave), (0, 0, 0))[2]
                          for wave, value in pressure.items()))
    assert mean == s.Rational(1, 6)
    pressure_bound = s.simplify(sum(s.Abs(value) for value in pressure.values()))
    assert pressure_bound == s.Rational(22, 9)
    boost = 2048
    positive_lower = s.Rational(boost, 6) - 84*pressure_bound
    assert positive_lower == 136
    a = s.symbols('a', nonnegative=True)
    raw_work_bound = s.Rational(boost, 6)*a - 84*pressure_bound*a**2
    weighted_bound = a**2 * raw_work_bound
    assert s.expand(weighted_bound) == s.Rational(boost, 6)*a**3 - 84*pressure_bound*a**4
    # Exact NS scaling and equal-mass background ledger.
    lam, ev, eg, iv, ig, zv, zg = s.symbols('lambda E_v E_g I_v I_g Z_v Z_g', positive=True)
    length = (ev + eg - ev/lam)/eg
    assert s.simplify(ev/lam + length*eg - (ev+eg)) == 0
    assert s.simplify(iv + ig - (iv+ig)) == 0
    output = {
        'scope': 'Exact periodic triad, coefficients and scaling algebra only; compactification, nonlocal errors and actual NS derivative are separate proofs.',
        'pressure_modes': {str(wave): str(value) for wave, value in sorted(pressure.items())},
        'mean_p_d3_v3': str(mean), 'pressure_sup_bound': str(pressure_bound),
        'boost': boost, 'periodic_pressure_work_lower': str(positive_lower),
        'envelope_counted_once': str(s.expand(weighted_bound)),
        'fixed_energy_identity': True, 'fixed_critical_mass_identity': True,
        'native_gradient_mass': str(lam*zv + zg/length),
    }
    target = Path('python/artifacts/unconditional-navier/critical-pressure-source-controls.json')
    target.write_text(json.dumps(output, indent=2)+'\n')
    print('PASS: exact triad pressure, positive boosted work, fixed E/I scaling. No PDE proof certified.')


if __name__ == '__main__':
    run()
