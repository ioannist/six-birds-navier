# Best-bound sheet for Cap(j) (attempt)

## (A) Setup (definitions)
- Shells S_j = {k: 2^j ≤ |k| < 2^(j+1)}, projections P_j, low-pass G_j = P_{≤j}.
- Leray projector: P_Leray; SGS closure:
  C_{≤j}(u) := G_j P_Leray ∇·(u_{≤j}⊗u_{≤j} - (u⊗u)_{≤j}).
- Port output: y_j := P_j C_{≤j}(u); port input: u_j := P_j u.
- Power: p_j(t) := ⟨u_j(t), y_j(t)⟩_{L^2}.
- Throughput target: ∫_s^t ⟨u_j, y_j⟩ dτ ≤ Λ(j;u) ∫_s^t |u_j|_2^2 dτ.

## (B) Attempt 1: scale-uniform bound with classical obstruction
- Commutator/para-product structure:
  C_{≤j} involves ∇·(u_{≤j}⊗u_{>j} + u_{>j}⊗u_{≤j}) + commutator terms.
- Apply product/commutator estimates to bound:
  ||P_j C_{≤j}(u)||_2 ≲ ||∇u||_∞ · ||u_j||_2.
- Obstruction term:
  F(t) := ||∇u(t)||_∞ (classical sup-norm obstruction channel; see obstruction_certificate.md).
- Final bound line (Attempt 1):
  |⟨u_j, y_j⟩| ≲ ||∇u||_∞ · ||u_j||_2^2.
- Energy inequality does not control F(t); this is the classical obstruction.

## (C) Attempt 2: energy-only bound (expect bad j-scaling)
- Use Bernstein on low-pass fields (3D):
  ||u_{≤j}||_∞ ≲ 2^(3j/2) ||u||_2,
  ||∇u_{≤j}||_∞ ≲ 2^(5j/2) ||u||_2.
- By the energy inequality, ||u(t)||_2 ≤ ||u(0)||_2, so the bound can be written in terms of ||u_0||_2.
- Substitute into the obstruction term to eliminate ||∇u||_∞:
  |⟨u_j, y_j⟩| ≲ 2^(5j/2) ||u||_2 · ||u_j||_2^2.
- Final j-scaling (Attempt 2):
  Λ(j) ≲ C · 2^(5j/2) · ||u||_2.
- DIV check: if Λ(j) grows like 2^(cj) with c>0, then ∑ 1/Λ(j) converges;
  energy-only bounds cannot yield No-Zeno.

## (D) Binary exit conclusion
- Outcome A (on track): Λ(j) ≲ (j+1)^α with α ≤ 1 from always-true estimates ⇒ No-Zeno feasible.
- Outcome B (obstruction found): Λ(j) depends on ||∇u||_∞ and energy-only fallback is exponential in j.
- Selected outcome: B (obstruction term is ||∇u||_∞; energy-only bound is exponential in j).
