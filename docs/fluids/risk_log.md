# Fluids risk log

## NS-T3 proof skeleton risks
- Mean-zero assumption required for inverting Δ / defining A^{-1}.
- Cancellation must use divergence-free field and Leray projection.
- Sobolev exponents must be correct (3D: ∇u in L∞ needs s>5/2).
- The previous high-norm-only interpolation produced a superlinear scalar ODE
  and did not close. The repaired proof uses the native-energy integrable
  coefficient ||Λ^{5/4}u||₂², finite-exponent periodic Sobolev embeddings and
  derivative product interpolation; μ must be strictly positive.
- Distinguish A=-PΔ from bare -Δ; domain restrictions matter.
- Inference metrics can be unstable with few shells (small N_co); avoid over-interpreting high-k fits.
- Steady/forced regime inference is still limited; larger averaging windows may be needed.
- Next step: strengthen Theorem-2 citation backbone and expand stationary inference coverage.

## Blow-up track risks
- We use positive supplied work W^+ to avoid cancellation.
- The frontier is defined only via storage growth since t_j.
- Throughput is an integrated inequality plus a feasibility budget; bounding Cap(j) growth is open.
- Route mismatch must be computed on packaged bridge objects, not restriction idempotence.
- The prior full-SGS shell-only strain estimate is false; the exact high–high
  counterexample and correct absolute power bound are in cap_bound_attempt.md.
- A loose exponential upper capacity bound fails to certify DIV but is not a
  lower bound on the minimal capacity or a no-go for all energy estimates.
- Anti-localization may fail in worst-case NS (intermittency); it must be stress-tested at higher N/Re and for adversarial initial data, and the Clay chain depends on proving HL-P2-ANTILOC structurally (not only observing it in forced runs).
