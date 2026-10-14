# (A) Objects and definitions
- Nonnegative effective capacity sequence: C(n) ≈ Cap(Pack_{0←n}) ≥ 0, n ∈ ℕ.
- Local step factor r(n) (single reduction multiplier).
- Route mismatch error e(n) (non-associativity penalty).

# (B) Approximate associativity ⇒ recursion
- Core recursion (route-stability abstraction):
  C(n+1) ≤ r(n)·C(n) + e(n).

# (C) Near-identity step model (tame reduction)
- Tame step bound:
  0 ≤ r(n) ≤ ((n+2)/(n+1))^p, with p ∈ ℕ.
- Telescoping product gives polynomial growth: ∏_{k< n} r(k) ≤ (n+1)^p.

# (D) Summable mismatch assumption
- Bounded partial sums:
  ∑_{k=0}^{n-1} e(k) ≤ E for all n.

# (E) Theorem (SBT hinge)
**Theorem (route-stability ⇒ polynomial capacity).**
If C(n) ≥ 0 for every n, e(n) ≥ 0, C(n+1) ≤ r(n)·C(n)+e(n),
r(n) ≤ ((n+2)/(n+1))^p, and ∑_{k<n} e(k) ≤ E, then
  C(n) ≤ (C(0)+E)·(n+1)^p.
In particular, capacity growth is not exponential in depth.

The nonnegative capacity assumption licenses replacing r(n)C(n) by the tame
upper factor times C(n). The Lean wrappers
`route_capacity_poly_bound_of_tame_factor` and
`route_capacity_poly_bound_of_tame_factor_of_bounded_sum` now include this
adapter explicitly. Nonnegative r is also needed for the stated product bound.

Mapping to NS route objects: fix shell j and set C(n) ↦ C(j,n), e(n) ↦ e(j,n); the factor r(n) ≤ ((n+2)/(n+1))^p corresponds to a tame per-depth inflation of Λ-certificates in the truncation ladder (route_objects_ns.md).

# (F) Plug-in to the No-Zeno engine
- If Cap(j) ≤ (j+1)^p with p ≤ 1, then ∑ 1/Cap(j) diverges ⇒ No-Zeno (constant w=θ).
- If p > 1, this upper bound alone does not decide DIV: the comparison series
  converges, but Cap(j) might be constant. Two-sided comparability
  Cap(j) ≍ (j+1)^p would give convergence. Positivity of Cap is required throughout.

# (G) What remains open (NS-specific)
- Prove route mismatch summability: ∑ e(n) ≤ E from packaging rules.
- Prove near-identity step factor r(n) ≤ ((n+2)/(n+1))^p under a fixed normalization.
- Assess whether r(n), e(n) reduce to ||∇u||_∞ control (classical obstruction risk).
