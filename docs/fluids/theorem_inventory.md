# Fluids theorem inventory

## Status table
| ID | Claim | Python evidence | Lean evidence | Status |
|---|---|---|---|---|
| NS-T1 | Leading closure is Laplacian-like (IR) | python/scripts/ns09_invariance_sanity.py; python/scripts/ns08_infer_closure_kernel.py; python/scripts/ns20_infer_closure_kernel_3d.py; python/scripts/ns21_infer_closure_kernel_3d_forced_avg.py; python/scripts/ns22_markov_lens_endomap_demo.py (Track B toy); artifacts: ns08_closure_kernel.png, ns09_invariance.png, ns20_closure_kernel_3d.png, ns21_closure_kernel_3d_forced_avg.png, ns22_markov_laplacian_scaling.png | formal/Fluids/Multiplier.lean (passivity ⇒ ν ≥ 0) | PARTIAL |
| NS-T2 | Energy inequality / accounting (abstract core + PDE ledger) | python/scripts/ns05_energy_ledger.py (ns05_energy_ledger.png); python/src/nswave/galerkin_core.py (demo) | formal/Fluids/Energy.lean (inner_rhs_self, inner_rhs_self_le_zero) | DONE |
| NS-T3 | Global regularity of completed layer (hyperviscous completion) | python/scripts/ns06_uv_completion_turns_on.py (ns06_uv_completion.png); python/scripts/ns14_energy_ledger_3d.py (ns14_energy_ledger_3d.png); python/scripts/ns15_uv_completion_turns_on_3d.py (ns15_uv_completion_3d.png); python/scripts/ns18_exponent_sanity_harness.py (ns18_exponent_sanity.png); proof dossier: docs/fluids/proof_skeleton_ns_alpha.md | none | PARTIAL (Proof skeleton complete + numeric stress tests; full write-up pending citations) |
| NS-T4 | μ→0 recovery on resolved scales | python/scripts/ns07_mu_to_zero_convergence.py (ns07_mu_convergence.png); python/scripts/ns16_mu_to_zero_convergence_3d.py (ns16_mu_convergence_3d.png) | formal/Fluids/Difference.lean | DONE (numeric) |
| NS-HOOK-B | Micro→macro Markov substrate + lens endomap shows low-k λ(k) ~ k^2 | python/scripts/ns22_markov_lens_endomap_demo.py (ns22_markov_laplacian_scaling.png) | none | HOOK / OPTIONAL |

- cite Lions / classical hyperviscous NS result (fill later).

## Open gaps
- 3D closure inference exists (NS-20/21), but high-k characterization is still noisy at small N_co; stationary turbulence regime not yet explored.
- Forcing/long averaging improves stability (NS-21), but broader parameter sweeps and stationary regimes remain.
- High-k steepening metrics are diagnostic-only at small N_co; avoid claiming a fitted p.
- Stronger mathematical bridge from empirical ℓ(k) to theorem statement not mechanized.
