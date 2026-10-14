#!/usr/bin/env python3
"""Exact rational finite-mode controls; never evidence for an NS trajectory."""
import argparse
from fractions import Fraction
import json
from pathlib import Path


def controls():
    gamma, cutoff = 3, 4
    rows = []
    for m in [1, 4, 5, 8, 16, 32]:
        # Conjugate transverse modes at +/-(m,0,0), amplitudes (0,1,0).
        kinetic = Fraction(2)
        strong = kinetic * (1 + m * m) ** gamma
        low = strong if m <= cutoff else Fraction(0)
        high = strong - low
        assert strong == low + high
        tail_fraction = high / strong
        if m > cutoff:
            assert tail_fraction == 1
        assert low <= kinetic * (1 + cutoff * cutoff) ** gamma
        weighted = strong * (1 + m * m)
        assert high <= weighted / (1 + cutoff * cutoff)
        rows.append({"frequency": m, "kinetic_energy": str(kinetic),
                     "strong_energy_gamma3": str(strong), "cutoff": cutoff,
                     "tail_fraction": str(tail_fraction),
                     "pointwise_cutoff_m_tail_fraction": "0"})
    low, high, q = Fraction(2), Fraction(3), Fraction(3, 5)
    total = low + high
    assert high <= q * total and total <= low / (1 - q)
    return {"status": "exact finite controls passed", "domain": "finite real transverse Fourier modes on a periodic control carrier",
            "original_target": "whole-space actual Navier maximal trajectory: OPEN",
            "scope": "These controls illustrate energy/strong-norm separation and cutoff drift. They are not NS trajectories or whole-space Clay instances.",
            "finite_modes": rows, "absorption_equality": str(total),
            "energy_class_spacetime_exponent": "10/3 (requires classical Sobolev/interpolation law)",
            "candidate_critical_exponent": "5 (independent supplier OPEN)"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = controls()
    text = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    print(result["status"])


if __name__ == "__main__":
    main()
