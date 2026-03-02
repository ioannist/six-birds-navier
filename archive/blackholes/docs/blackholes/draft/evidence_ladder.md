# BH Evidence Ladder (draft)

This section governs what we are allowed to claim in the BH write-up.
We separate algebraic theorems, numerical validation, time-domain realizability, and hooks.
Claims should follow the tier boundaries below.

## Tier 1 — Algebraic / theorem-level (Lean-backed where available)
- What we can claim:
  - BH-T2 passivity => |R|<=1 via Cayley transform; Lean-backed.
- What we cannot claim yet:
  - Flux identity j=omega(|R|^2-1) is not mechanized in Lean.
- Evidence pointers:
  - formal/BlackHoles/Impedance.lean (lemma abs_R_from_Z_le_one_of_re_nonneg)
  - python/tests/test_bh04_passivity_equivalence.py
  - python/scripts/bh04_passivity_checks.py

## Tier 2 — Frequency-domain numerical validation (ODE scattering)
- What we can claim:
  - Scattering solver matches analytic free case; echo formula matches full ODE solver.
  - Frequency-dependent Debye Z(omega) validated; RW barrier computed and compared to PT.
- What we cannot claim yet:
  - No rigorous convergence proof for ODE solver beyond numeric agreement.
- Evidence pointers:
  - python/tests/test_bh03_free_case.py
  - python/scripts/bh03_free_case_demo.py
  - python/scripts/bh07_echo_formula_validation.py
  - python/tests/test_bh07_echo_formula_validation.py
  - python/scripts/bh15_echo_freqdep_validation.py
  - python/tests/test_bh15_echo_freqdep_validation.py
  - python/scripts/bh17_rw_barrier_compare.py
  - python/artifacts/bh17_rw_vs_pt.png

## Tier 3 — Time-domain realizability + stability/energy ledger
- What we can claim:
  - Realizable boundary via auxiliary states produces echo waveform.
  - Discrete energy ledger for passive boundary is non-injecting within tolerance.
- What we cannot claim yet:
  - Full stability region for dt/dx and parameter ranges is not mapped.
- Evidence pointers:
  - python/scripts/bh08_time_domain_demo.py
  - python/artifacts/bh08_waveform.png
  - python/scripts/bh13_energy_ledger_demo.py
  - python/artifacts/bh13_energy.png
  - python/scripts/bh14_boundary_stability_compare.py
  - python/artifacts/bh14_energy_compare.png

## Tier 4 — Hooks / out-of-scope regimes (explicitly labeled)
- What we can claim:
  - Superradiance ledger is a toy hook; gain > 1 leads to instability in the toy model.
- What we cannot claim yet:
  - Kerr/Teukolsky physics not implemented; only toy active-band reflectivity.
  - Full positive-real/analyticity mechanization in Lean is not attempted; KK is numeric.
- Evidence pointers:
  - python/scripts/bh11_superradiance_ledger_demo.py
  - python/artifacts/bh11_superradiance_gain.png
  - docs/blackholes/risk_log.md
  - python/scripts/bh05_kramers_kronig_numeric.py
  - python/scripts/bh10_cross_channel_linkage.py

## Paste-ready paragraph outline
- Tier 1 (P6 accounting): Lean-backed Cayley inequality underwrites the non-superradiant ledger claim; we cite only the mechanized algebra.
- Tier 2 (P5 packaging): frequency-domain solver evidence validates the packaged boundary response and the echo transfer formula (including dispersive Z).
- Tier 3 (P6 accounting + realizability): time-domain auxiliary-state boundary produces echoes and respects the discrete energy ledger within tolerance.
- Tier 4 (hooks): superradiance is treated as a toy ledger extension; Kerr/Teukolsky and full positive-real Lean formalization remain out of scope.
