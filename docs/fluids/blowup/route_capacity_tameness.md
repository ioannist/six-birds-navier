# (A) Objects and definitions
- Effective capacity sequence: C(n) ≈ Cap(Pack_{0←n}), n ∈ ℕ.
- Local step factor r(n) (single reduction multiplier).
- Route mismatch error e(n) (non-associativity penalty).

# (B) Approximate associativity ⇒ recursion
- Core recursion (route-stability abstraction):
  C(n+1) ≤ r(n)·C(n) + e(n).

# (C) Near-identity step model (tame reduction)
- Tame step bound:
  r(n) ≤ ((n+2)/(n+1))^p, with p ∈ ℕ.
- Telescoping product gives polynomial growth: ∏_{k< n} r(k) ≤ (n+1)^p.

# (D) Summable mismatch assumption
- Bounded partial sums:
  ∑_{k=0}^{n-1} e(k) ≤ E for all n.

# (E) Theorem (SBT hinge)
**Theorem (route-stability ⇒ polynomial capacity).**
If C(0) ≥ 0, e(n) ≥ 0, C(n+1) ≤ r(n)·C(n)+e(n),
r(n) ≤ ((n+2)/(n+1))^p, and ∑_{k<n} e(k) ≤ E, then
  C(n) ≤ (C(0)+E)·(n+1)^p.
In particular, capacity growth is not exponential in depth.

Mapping to NS route objects: fix shell j and set C(n) ↦ C(j,n), e(n) ↦ e(j,n); the factor r(n) ≤ ((n+2)/(n+1))^p corresponds to a tame per-depth inflation of Λ-certificates in the truncation ladder (route_objects_ns.md).

# (F) Plug-in to the No-Zeno engine
- If Cap(j) ≤ (j+1)^p with p ≤ 1, then ∑ 1/Cap(j) diverges ⇒ No-Zeno (constant w=θ).
- If p > 1, DIV fails even with polynomial capacity; hinge is “prove p ≤ 1 or improve bound.”

# (G) What remains open (NS-specific)
- Prove route mismatch summability: ∑ e(n) ≤ E from packaging rules.
- Prove near-identity step factor r(n) ≤ ((n+2)/(n+1))^p under a fixed normalization.
- Assess whether r(n), e(n) reduce to ||∇u||_∞ control (classical obstruction risk).
