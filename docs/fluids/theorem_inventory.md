# Fluids theorem inventory

## Status table
| ID | Claim | Python evidence | Lean evidence | Status |
|---|---|---|---|---|
| NS-T1 | Scalar time-local IR implication under explicit C² and right-neighborhood passivity premises | docs/fluids/passive_low_frequency_closure.md (analytic derivative/Taylor return); python/scripts/ns09_invariance_sanity.py; python/scripts/ns08_infer_closure_kernel.py; python/scripts/ns20_infer_closure_kernel_3d.py; python/scripts/ns21_infer_closure_kernel_3d_forced_avg.py; python/scripts/ns22_markov_lens_endomap_demo.py (Track B toy) | formal/Fluids/Multiplier.lean (right-neighborhood passivity ⇒ η ≥ 0; fixed-mode counterexample) | Analytic conditional scalar-symbol implication supplied; constitutive applicability remains PARTIAL |
| NS-T2 | Energy cancellation and dissipative inner-product algebra | python/scripts/ns05_energy_ledger.py (ns05_energy_ledger.png); python/src/nswave/galerkin_core.py (demo) | formal/Fluids/Energy.lean (inner_rhs_self, inner_rhs_self_le_zero) | DONE (abstract algebra); continuous PDE energy balance not mechanized |
| NS-T3 | Global smoothness of hyperviscous completion for fixed μ>0, ν≥0, α≥5/4 | python/scripts/ns06_uv_completion_turns_on.py; python/scripts/ns14_energy_ledger_3d.py; python/scripts/ns15_uv_completion_turns_on_3d.py; python/scripts/ns18_exponent_sanity_harness.py; corrected analytic proof: docs/fluids/proof_skeleton_ns_alpha.md | none | Analytic critical-energy proof supplied; not mechanized, no μ-uniform estimate |
| NS-T4 | Conditional μ→0 recovery with common data/force on a supplied smooth original NS interval | docs/fluids/smooth_interval_recovery.md (analytic O(μ) proof with explicit derivative cost and bootstrap); python/scripts/ns07_mu_to_zero_convergence.py; python/scripts/ns16_mu_to_zero_convergence_3d.py | formal/Fluids/Difference.lean (abstract derivative/forcing algebra only) | Analytic conditional theorem restored; no unconditional original NS continuation or Lean PDE limit theorem |
| NS-HOOK-B | Micro→macro Markov substrate + lens endomap shows low-k λ(k) ~ k^2 | python/scripts/ns22_markov_lens_endomap_demo.py (ns22_markov_laplacian_scaling.png) | none | HOOK / OPTIONAL |

- The repaired hyperviscous proof gives the explicit periodic/forcing adapters and primary source link.

## Open gaps
- 3D closure inference exists (NS-20/21), but high-k characterization is still noisy at small N_co; stationary turbulence regime not yet explored.
- Forcing/long averaging (NS-21) is an empirical diagnostic; improved fits or stability are not guaranteed. Empty retained frequency bands must remain unavailable.
- High-k steepening metrics are diagnostic-only at small N_co; avoid claiming a fitted p.
- Stronger mathematical bridge from empirical ℓ(k) to theorem statement not mechanized.
