# Divergence decision (No-Zeno engine)

## A. Decision engine (constant work quantum)
- Work quantum per crossing: W_j^+[t_j,t_{j+1}] ≥ w(j), with w(j)=θ (constant).
- Throughput/budget: W_j^+[t_j,t_{j+1}] ≤ Cap(j)·Δt_j.
- Hence Δt_j ≥ w(j)/Cap(j) = θ/Cap(j).
- If partial sums ∑_{k< n} w(k)/Cap(k) are unbounded, then t_n → ∞ (Lean: Fluids.ZenoWorkCap).
- With constant w=θ, the divergence target is ∑ 1/Cap(j) = ∞ (up to constants).
- Cap(j) here is the capacity that upper-bounds positive work W_j^+.

## B. Divergence facts (checkable)
1) Exponential growth:
   - If Cap(j) ≳ 2^{β j} (β>0), then ∑ 1/Cap(j) converges (geometric series).
   - No-Zeno cannot be proved with constant w if Cap grows exponentially.
2) Polynomial growth:
   - If Cap(j) ≲ (j+1)^α, then ∑ 1/Cap(j) diverges iff α ≤ 1.
3) Time-dependent obstruction factor:
   - If Cap(j) ≲ F(t)·G(j), divergence requires control on F over [t_j,t_{j+1}].
   - With F(t)=||∇u||_∞ (or ||ω||_∞), this is exactly the classical obstruction channel.

## C. Apply to NS-BU-14 bounds
- Attempt 1 (best commutator bound):
  |⟨u_j, y_j⟩| ≲ ||∇u||_∞ ||u_j||_2^2, so Cap(j) depends on F(t)=||∇u||_∞.
  DIV cannot be decided without controlling the classical obstruction norm.
- Attempt 2 (energy-only + Bernstein):
  Λ(j) ≲ C·2^{5j/2}·||u||_2, hence Cap(j) grows exponentially in j.
  Therefore ∑ 1/Cap(j) converges and DIV fails for constant w=θ.

## D. Decision / next lemma slots
- Outcome: DIV is not certifiable from energy inequality + LP calculus alone.
- Next lemma options:
  - HL-CAP upgrade: prove Cap(j) grows at most polynomial (≤ linear) without ||∇u||_∞.
  - HL-WORK upgrade: allow w(j) to grow with j (nontrivial to justify).
  - Scale-passivity upgrade: derive bounded dissipation density / bounded SGS response by shell.
