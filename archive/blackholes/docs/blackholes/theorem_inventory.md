### Status table
| ID | Claim (short) | Python evidence | Lean evidence | Status |
|---|---|---|---|---|
| BH-T1 | Boundary response completeness (linear exterior) | BH-03 solver (parameterized by Z/R), BH-02 conventions sanity | none | PARTIAL |
| BH-T2 | Passivity ⇒ |R|≤1 (non-superradiant) | bh04_passivity_checks.py + tests | abs_R_from_Z_le_one_of_re_nonneg | DONE |
| BH-T3 | Echo/cavity transfer function | bh07_echo_formula_validation.py, bh08_time_domain_demo.py | none | DONE (numeric) |
| BH-T4 | Dispersion links conservative/dissipative | bh05_kramers_kronig_numeric.py, bh10_cross_channel_linkage.py | none | DONE (numeric) |
| BH-HOOK | Superradiance ledger modifies inequality | bh11_superradiance_ledger_demo.py | none | HOOK |

## BH-T1 Boundary response completeness
- Statement (reduced 1D model, fixed real ω>0):
  For the exterior ODE ψ'' + (ω^2 - V(x))ψ=0 on x∈[x0,∞), the space of local linear boundary conditions at x=x0 is equivalent to specifying the complex ratio
    Y(ω) := ψ'(x0,ω) / ψ(x0,ω)
  (or normalized impedance Z(ω) := iY(ω)/ω, or reflectivity R(ω) := A_out/A_in).
- Assumptions:
  (i) linear time-invariant interior coupling as seen from the exterior,
  (ii) local boundary relation at x0 in frequency domain (no explicit dependence on distant future),
  (iii) ω chosen away from singular points where ψ(x0)=0 for the selected solution.
- Evidence / implementation:
  BH-02 fixes a consistent convention and maps (Z↔R↔Y).
  BH-03 implements scattering_reflection(ω,V,x0,x_max, Z or R_boundary) and matches the V≡0 analytic case at <1e-6.
- Mechanization:
  Not yet in Lean (we currently mechanize the core Cayley inequality used in BH-T2).

## BH-T2 Passivity ⇒ |R|≤1
- Statement (scalar 1D, real omega > 0, non-superradiant):
  passive inner boundary (no net energy injection into exterior) implies |R(omega)| <= 1.
- Conventions: j := Im(psi* psi') (flux in +x).
  For psi = e^{-i omega y} + R e^{+i omega y}, j = omega (|R|^2 - 1).
  Thus passivity (j(0) <= 0) <=> |R| <= 1.
- Boundary mapping: R = (1 - Z) / (1 + Z).
  For real omega, |R| <= 1 <=> Re(Z) >= 0 via
  |R|^2 = |1 - Z|^2 / |1 + Z|^2 and |1 + Z|^2 - |1 - Z|^2 = 4 Re(Z).
- Mechanization hook: Cayley-transform inequality in Complex;
  see formal/BlackHoles/Impedance.lean lemma.
- Evidence: python/scripts/bh04_passivity_checks.py, python/tests/test_bh04_passivity_equivalence.py
- Lean: formal/BlackHoles/Impedance.lean lemma abs_R_from_Z_le_one_of_re_nonneg

## BH-T3 Echo / cavity transfer function
- Cavity transfer factor: H(omega) = 1 / (1 - R0(omega) * r(omega) * exp(2 i omega Delta x)).
- Predicted reflection: R_pred(omega) = r(omega) * exp(-2 i omega Delta x)
  + (t(omega)^2 * R0(omega)) / (1 - R0(omega) * r(omega) * exp(2 i omega Delta x)).
- Validated numerically against the full ODE scattering solver in
  python/scripts/bh07_echo_formula_validation.py for |R0| <= 0.2.
- Time-domain FDTD demo: python/scripts/bh08_time_domain_demo.py produces an echo train with
  spacing ≈ 2*Delta x (Delta x = x_ref - x0).
- Evidence: bh07_echo_formula_validation.py matches full ODE solver at ~1e-7 abs; bh08_time_domain_demo.py shows echo spacing ≈2Δx.
- Also validated for frequency-dependent Debye Z(ω): python/scripts/bh15_echo_freqdep_validation.py.
- End-to-end FDTD→fit demo: python/scripts/bh16_fdtd_inference_demo.py recovers Debye Z(ω) using a two-point plane-wave estimator.

## BH-T4 Cross-channel dispersion linkage
- Conservative proxy: Re Z(0) (Love-like static response in the toy model).
- Dissipative proxy: Im Z(omega).
- KK link at omega -> 0:
  Re Z(0) = re_inf + (2/pi) ∫_0^∞ Im Z(omega) / omega d omega.
- Validated numerically in python/scripts/bh10_cross_channel_linkage.py.
- Evidence: BH-05 numeric KK reconstruction; BH-10 static KK linkage ReZ(0) from ImZ/ω, errors ~1e-4.

### Open gaps
- Kerr/Teukolsky (mode mixing, superradiance) not implemented; BH-11 is a toy ledger hook only.
- Full causal/positive-real formalization in Lean not attempted (numeric-only KK checks).
- Time-domain impedance boundary stability: dt reduced (CFL 0.6) for long runs; investigate stronger stability criterion if needed.
