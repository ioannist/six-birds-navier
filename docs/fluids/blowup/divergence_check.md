# Divergence conditions for the No Zeno engine

Assume θ>0, Cap(j)>0 and θ≤Cap(j)Δt_j on actual crossing intervals.
Then Δt_j≥θ/Cap(j). Unbounded partial sums of θ/Cap(j) exclude a finite
upper bound on all crossing times (`Fluids.ZenoWorkCap`).

## Correct comparison directions

- Cap(j)≤K(j+1)^α eventually, K>0, α≤1: the reciprocal series diverges.
- Cap(j)≥k(j+1)^α eventually, k>0, α>1: the reciprocal series converges.
- Cap(j)≍(j+1)^α: divergence holds exactly when α≤1.
- Cap(j)≥k 2^{βj}, k,β>0: the reciprocal series converges.
- An exponential upper bound alone decides nothing about the actual series.
  The counterexample Cap(j)=1 satisfies such an upper bound and has divergent
  reciprocal sum. The same counterexample invalidates an “iff” conclusion
  from only a polynomial upper bound with α>1.

A convergent lower bound on latency merely leaves No-Zeno unproved. It does
not show actual crossing increments have finite sum or that blow-up occurs.
For a time-dependent power majorant F(t), one must control its integral on
the actual intervals; time integrability is not a uniform pointwise bound.

## Actual SGS bound

[cap_bound_attempt.md](cap_bound_attempt.md) repairs the full SGS estimate:
energy gives an absolute bound |p_j|≤C2^{5j/2}||u0||₂³ for unforced NS,
not Λ(j)≤C2^{5j/2}||u0||₂ multiplying shell input energy. Choosing that
loose absolute Cap(j) gives a convergent reciprocal comparison. The minimal
capacity and its DIV status remain unresolved. No impossibility theorem for
all stronger source estimates follows from this failure.
