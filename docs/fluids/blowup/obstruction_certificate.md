# Scope of the discrete ICAP obstruction

For a finite scalar gain profile g, the hypothesis

    max(Σ_i g(i)u(i)²,0) ≤ Λ Σ_i u(i)²

on a feasible class containing every coordinate basis vector implies
max(g(i),0)≤Λ for every i, by testing that vector. This is the exact
mechanized no-go in `formal/Fluids/IcapNoGo.lean`; the restricted-feasibility
version now states basis admissibility explicitly. It needs no PDE estimate.

## What the nonlinear SGS port does not inherit

The actual nonlinear bridge and its supplied-force sign are fixed in
[ns_bridge_object.md](ns_bridge_object.md). The previously asserted bound
|p_j|≲||∇u||∞||u_j||₂² for its full SGS output is false. The exact high–high
triad and a correct absolute power bound are derived in
[cap_bound_attempt.md](cap_bound_attempt.md).

Consequently the discrete diagonal-gain theorem cannot simply be applied
to that nonlinear, history-dependent bridge: an actual gain representation,
solenoidal localized test feasibility, and a return to the claimed strain
readout would first be needed. Finite-scale localized averages also require
a limit and fixed-background control before they imply a pointwise bound.
These adapters remain OPEN. The discrete theorem by itself is not a
continuum ICAP-to-strain theorem or a no-go for all energy methods.

## Valid No Zeno return and remaining obligations

If actual crossings satisfy w(j)≤Cap(j)(t_{j+1}-t_j), Cap(j)>0, and
Σ_j w(j)/Cap(j) has unbounded partial sums, the Lean No-Zeno theorem applies.
The correct energy-only absolute CAP bound grows like 2^{5j/2} for fixed
native energy. Its reciprocal comparison does not certify DIV with constant
work. This does not prove that no sharper capacity could work.

Anti-localization can exclude the specific localized tests used by the
obstruction. It neither supplies ICAP nor forces WORK, DIV, persistence
of feasibility, or a continuation criterion. No-Zeno implies no blow-up
only after the explicitly open storage-frontier/PDE return is proved.
The conditional original scaffold is retained with these obligations visible.
