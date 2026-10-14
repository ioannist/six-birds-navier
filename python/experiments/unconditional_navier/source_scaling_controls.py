#!/usr/bin/env python3
"""Exact scaling diagnostics for profiles, not Navier solution trajectories."""
from fractions import Fraction
import json
from pathlib import Path


def diagnostic(alpha, beta):
    """u = tau**(-beta) U(x/tau**alpha), tau approaches zero."""
    powers = {
        "enstrophy": 2 * beta - alpha,
        "squared_enstrophy_dose": 4 * beta - 2 * alpha,
        "L5_dose": 5 * beta - 3 * alpha,
        "Wiener_squared_dose": 2 * beta,
        "positive_middle_strain_dose": beta + alpha,
    }
    return {
        "alpha": str(alpha), "beta": str(beta),
        "kinetic_energy_decay_power": str(3 * alpha - 2 * beta),
        "bounded_kinetic_energy": 3 * alpha - 2 * beta >= 0,
        "doses": {name: {"blowup_power": str(power),
                         "time_integrable": power < 1,
                         "logarithmic_borderline": power == 1}
                  for name, power in powers.items()},
    }


def main():
    cases = {
        "parabolic_selfsimilar_profile": diagnostic(Fraction(1, 2), Fraction(1, 2)),
        "historical_energy_profile": diagnostic(Fraction(4, 9), Fraction(2, 3)),
    }
    parabolic = cases["parabolic_selfsimilar_profile"]
    assert parabolic["bounded_kinetic_energy"]
    assert parabolic["doses"]["enstrophy"]["time_integrable"]
    assert all(d["logarithmic_borderline"] for k, d in parabolic["doses"].items()
               if k != "enstrophy")
    historical = cases["historical_energy_profile"]
    assert historical["kinetic_energy_decay_power"] == "0"
    assert historical["doses"]["enstrophy"]["time_integrable"]
    assert not any(d["time_integrable"] for k, d in historical["doses"].items()
                   if k != "enstrophy")
    result = {"scope": "Fixed nonzero Schwartz profile with nonzero positive-middle-strain readout; no assertion that a scaled profile solves the original unforced PDE.",
              "integrability_rule": "Integral of tau^(-p) near zero is finite exactly when p<1; p=1 diverges logarithmically.",
              "cases": cases}
    root = Path(__file__).resolve().parents[3]
    output = root / "python/artifacts/unconditional-navier/source-scaling-controls.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print("Exact scaling controls passed; profile diagnostics only:", output)


if __name__ == "__main__":
    main()
