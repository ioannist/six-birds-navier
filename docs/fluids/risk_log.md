# Fluids risk log

## NS-T3 proof skeleton risks
- Mean-zero assumption required for inverting Δ / defining A^{-1}.
- Cancellation must use divergence-free field and Leray projection.
- Sobolev exponents must be correct (3D: ∇u in L∞ needs s>5/2).
- Interpolation step requires α>1; threshold α≥5/4 used in absorption.
- Commutator estimate must be justified on T^3 (periodic setting).
- Distinguish A=-PΔ from bare -Δ; domain restrictions matter.
- Inference metrics can be unstable with few shells (small N_co); avoid over-interpreting high-k fits.
- Steady/forced regime inference is still limited; larger averaging windows may be needed.
- Next step: strengthen Theorem-2 citation backbone and expand stationary inference coverage.

## Blow-up track risks
- We use positive supplied work W^+ to avoid cancellation.
- The frontier is defined only via storage growth since t_j.
- Throughput is an integrated inequality plus a feasibility budget; bounding Cap(j) growth is open.
- Route mismatch must be computed on packaged bridge objects, not restriction idempotence.
- BU-16: canonical work-quantum scaling appears to require ||∇u||_∞ control (BKM-type channel).
- BU-16: without a new scale-passivity lemma, No-Zeno cannot be certified from energy-only bounds.
- Anti-localization may fail in worst-case NS (intermittency); it must be stress-tested at higher N/Re and for adversarial initial data, and the Clay chain depends on proving HL-P2-ANTILOC structurally (not only observing it in forced runs).
