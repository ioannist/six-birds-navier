# (A) Setup (canonical objects)
- Shells S_j = {k: 2^j ≤ |k| < 2^(j+1)}, projections P_j, low-pass G_j = P_{≤j}.
- Port input u_j := P_j u, port output y_j := P_j C_{≤j}(u).
- Power p_j := ⟨u_j, y_j⟩_{L^2}, positive work W_j^+[s,t] := ∫_s^t max(p_j,0) dτ.
- Capacity inequality: W_j^+[s,t] ≤ Λ(j) ∫_s^t |u_j|_2^2 dτ and Cap(j) := Λ(j) \bar B(j).

# (B) What “canonical w(j)” means in NS terms
- w(j) must be tied to NS dyadic objects (flux/structure), not arbitrary thresholds.
- Frontier definition should satisfy: blow-up (if any) implies j_*(t) → ∞ under smoothness/LP assumptions.
- w(j) should arise from a structural increase that forces positive work by passivity.

# (C) Candidate activation definitions and induced w(j)
1) Flux-activation candidate (interscale transfer):
   - Activation: A_{j+1}(t) := sup_{τ∈[t_j,t]} ∫_{t_j}^τ Π_j^+(s) ds,
     where Π_j is a dyadic flux proxy (symbolic).
   - Crossing: t_{j+1} := inf{t ≥ t_j : A_{j+1}(t) ≥ θ_j}.
   - Work quantum target: w(j) := θ_j with W_j^+[t_j,t_{j+1}] ≥ w(j).
2) Weighted-regularity candidate (structure growth):
   - Activation: A_{j+1}(t) := sup_{τ∈[t_j,t]} 2^{2s(j+1)} ||u_{j+1}(τ)||_2^2.
   - Crossing: t_{j+1} := inf{t ≥ t_j : A_{j+1}(t) ≥ Θ_j}.
   - Work quantum target: w(j) := Θ_j (structural threshold tied to H^s scale).

# (D) Attempted lower bound W_j^+ ≥ w(j)
- Flux candidate: need Π_j^+ to control ⟨u_j, y_j⟩_+; missing lemma:
  “positive interscale transfer forces positive SGS work at the port.”
- Weighted-regularity candidate: need a storage inequality linking
  2^{2s(j+1)}||u_{j+1}||_2^2 increase to W_j^+; missing lemma:
  “structure growth implies positive work on shell j.”
- Failure point: no canonical inequality from NS that lower-bounds W_j^+ by these activations without extra assumptions.

# (E) Scaling arithmetic vs Cap(j)
- Best known Cap growth (energy-only): Λ(j) ≲ C·2^{5j/2}·||u_0||_2, so Cap(j) grows exponentially.
- If w(j) is at most polynomial in j (or constant), then ∑ w(j)/Cap(j) converges (DIV fails).
- To cancel exponential Cap, w(j) would need exponential growth in j.
- Proving exponential w(j) would require a strong lower bound on interscale transfer,
  which appears to depend on ||∇u||_∞ (classical obstruction channel).

# (F) Obstruction certificate
- Certificate: any w(j) large enough to offset the best available Cap(j) bound
  requires controlling ||∇u||_∞ (or equivalently a BKM-type channel) to force
  sustained positive interscale transfer; this is beyond energy-only estimates.

# (G) Update hooks (lemma slots)
- NZ-L2/NZ-L3: need a proven link from activation to W_j^+.
- NZ-L4/NZ-L5: need a Cap(j) growth bound ≤ polynomial without ||∇u||_∞.
- BF-L2/BF-L6: identify which SP-* (scale-passivity) axiom fails first.
