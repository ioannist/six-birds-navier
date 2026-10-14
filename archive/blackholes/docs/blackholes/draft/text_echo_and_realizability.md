# Echo transfer and realizability (BH-T3)

## BH-T3: cavity transfer function (frequency domain)
- TODO(BH): state formula for R_out in terms of r, t, R0, Delta x; define Delta x = x_ref - x0
- TODO(BH): note assumptions (1D barrier, V->0 at x_max, symmetric barrier)

## Validation vs full ODE solver
- TODO(BH): cite BH-07 (constant R0) + BH-15 (freq-dependent Debye Z)
- TODO(BH): mention typical error scale (~1e-7) and where it is measured

## Time-domain realizability: auxiliary-state impedance boundary
- TODO(BH): write the two equations (psi_x = a0 psi_t + sum a_k q_k; q' = -b q + b psi_t)
- TODO(BH): connect to realizable passive boundary response

## Discrete energy ledger (P6 accounting certificate)
- TODO(BH): define discrete energy E including boundary storage term 0.5 sum (a_k/b_k) q_k^2
- TODO(BH): cite BH-13 demo + max_rel_increase threshold

## Implementation stability notes (non-theorem)
- TODO(BH): exact q-update, predictor/corrector, dt choice in BH-08, BH-14
