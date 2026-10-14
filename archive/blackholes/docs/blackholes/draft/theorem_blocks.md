Skeleton theorem blocks for BH section; each block lists claim, assumptions,
what is proved vs checked, and evidence pointers. Draft-only, no final prose.

## BH-T1 Boundary response completeness (Packaging / P5)
- Claim: Exterior ODE closes given Y(omega) or Z(omega) or R(omega) at x0.
  The parameterizations are equivalent and define the boundary response class.
- Assumptions:
  - Linear exterior ODE; V is fixed and real in the exterior.
  - Local boundary operator in the frequency domain at x0.
  - Omega chosen away from singularities where psi(x0,omega)=0.
- What is proved vs checked: algebraic equivalence of parameterizations; numeric solver sanity.
- Evidence pointers:
  - docs/blackholes/draft/notation.md
  - python/src/bhwave/conventions.py
  - python/src/bhwave/scattering.py
  - python/tests/test_bh03_free_case.py
  - python/scripts/bh03_free_case_demo.py

## BH-T2 Passivity => |R(omega)| <= 1 (non-superradiant ledger / P6)
- Claim: passivity implies |R|<=1; equivalently Re(Z)>=0.
- Assumptions: real omega>0; non-superradiant ledger; linear response.
- What is proved vs checked: Lean mechanizes Cayley inequality; flux identity stated; numeric checks.
- Evidence pointers:
  - docs/blackholes/theorem_inventory.md (BH-T2)
  - python/scripts/bh04_passivity_checks.py
  - python/tests/test_bh04_passivity_equivalence.py
  - formal/BlackHoles/Impedance.lean (lemma abs_R_from_Z_le_one_of_re_nonneg)

## BH-T3 Echo transfer function (cavity / Theorem-3 style)
- Claim: R_out = r e^{-2 i omega Dx} + (t^2 R0)/(1 - R0 r e^{2 i omega Dx}).
- Assumptions: 1D barrier + cavity; symmetric barrier; V->0 at x_max; free-wave matching.
- What is proved vs checked: algebraic derivation stated; validated vs ODE solver; validated for freq-dependent Z(omega).
- Evidence pointers:
  - python/src/bhwave/echo.py
  - python/scripts/bh07_echo_formula_validation.py
  - python/tests/test_bh07_echo_formula_validation.py
  - python/scripts/bh15_echo_freqdep_validation.py
  - python/tests/test_bh15_echo_freqdep_validation.py

## BH-T4 Dispersion links conservative/dissipative parts (KK / P6)
- Claim: ReZ reconstructed from ImZ via KK; static limit ReZ(0) = re_inf + (2/pi) * int ImZ/omega.
- Assumptions: causal analytic Z in upper half-plane; passive family used in demo.
- What is proved vs checked: numeric KK only (no Lean).
- Evidence pointers:
  - python/src/bhwave/kramers_kronig.py
  - python/scripts/bh05_kramers_kronig_numeric.py
  - python/tests/test_bh05_kramers_kronig.py
  - python/scripts/bh10_cross_channel_linkage.py
  - python/tests/test_bh10_cross_channel_linkage.py

## BH-HOOK Superradiance ledger and instability condition
- Claim: when gain |R0 r|>1, cavity pole has Im(omega)>0 (toy); non-superradiant |R|<=1 replaced by ledger-aware constraint.
- Assumptions: toy active-band R0(omega); Kerr details omitted.
- What is proved vs checked: toy algebra + demo plot.
- Evidence pointers:
  - python/src/bhwave/superradiance.py
  - python/scripts/bh11_superradiance_ledger_demo.py
  - python/tests/test_bh11_superradiance_algebra.py
  - docs/blackholes/risk_log.md
