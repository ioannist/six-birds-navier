# Target statement (Theorem 2)
- For α ≥ 5/4 on T^3, smooth divergence-free mean-zero u0 and smooth f give a unique global smooth solution.
- Periodic setting with Leray projector P and Stokes operator A = -PΔ.

# Setup / objects
- H := divergence-free, mean-zero L^2(T^3) vector fields.
- P Leray projector; A = -PΔ; B(u,v)=P((u·∇)v).
- P commutes with Λ^s (Fourier multipliers) on T^3.

# Scaling / criticality checkpoint
- Scaling: x→λx, t→λ^{2α}t, u→λ^{2α-1}u.
- In d=3, ||u_λ||_2 scales like λ^{2α-5/2}.
- α=5/4 is energy-critical; α≥(d+2)/4 gives subcritical control (d=3 ⇒ α≥5/4).

# Lemma chain inventory (names + templates)
1) L0: Nonlinear cancellation ⟨B(u,u),u⟩=0.
   - Where used: base energy identity; sets up dissipation ledger.
   - Reference stub: standard NS identity.
2) L1: Local existence/uniqueness in H^s (s>5/2).
   - Where used: initial local well-posedness + continuation.
   - Reference stub: Temam; semilinear evolution on T^3.
3) L2: Fractional Leibniz / Kato–Ponce commutator.
   - Where used: bound ⟨Λ^s(u·∇u),Λ^s u⟩.
   - Reference stub: Kato–Ponce.
4) L3: Sobolev embedding H^s→L∞ (s>3/2) and ∇u in L∞ for s>5/2.
   - Where used: control ||∇u||_∞.
   - Reference stub: standard Sobolev on T^3.
5) L4: Interpolation between H^s and H^{s+α}.
   - Where used: relate ||Λ^{s+1}u||_2 to ||Λ^{s+α}u||_2.
   - Reference stub: Gagliardo–Nirenberg.
6) L5: Young / absorption inequality.
   - Where used: absorb nonlinear term into μ||Λ^{s+α}u||_2^2.
   - Reference stub: standard energy method.
7) L6: Continuation criterion (e.g., ∫||∇u||_∞ finite ⇒ continuation).
   - Where used: global extension.
   - Reference stub: NS blow-up criterion.
8) L7: Grönwall / bootstrap closure.
   - Where used: close differential inequality on ||Λ^s u||_2.
   - Reference stub: standard ODE inequality.

# Energy-method skeleton (bullet proof flow)
- Apply Λ^s to the equation, take L^2 inner product with Λ^s u.
- Use L0 to cancel ⟨Λ^s B(u,u), Λ^s u⟩ at s=0.
- Write: (1/2)d/dt||Λ^s u||_2^2 + ν||Λ^{s+1}u||_2^2 + μ||Λ^{s+α}u||_2^2 ≤ |⟨Λ^s B(u,u), Λ^s u⟩| + forcing work.
- Bound nonlinear term via commutator: |⟨Λ^s(u·∇u), Λ^s u⟩| ≤ C ||∇u||_∞ ||Λ^s u||_2^2.
- Alternative bound: |⟨Λ^s(u·∇u), Λ^s u⟩| ≤ C ||Λ^{s+1}u||_2 ||Λ^s u||_2^2.
- Use L3 to control ||∇u||_∞ from H^s (s>5/2).
- Interpolate: ||Λ^{s+1}u||_2 ≤ ||Λ^s u||_2^{1-1/α} ||Λ^{s+α}u||_2^{1/α} (α>1).
- Apply Young to absorb a portion of ||Λ^{s+α}u||_2^2 into the μ term.
- THIS is where α≥5/4 ensures the exponent arithmetic yields a subcritical differential inequality.
- Obtain a closed differential inequality for ||Λ^s u||_2^2 with forcing term.
- Use Grönwall to bound ||Λ^s u||_2 uniformly on any finite interval.
- Invoke L6 continuation criterion to extend globally.

# Proof stress-test checklist (things that can go wrong)
- Missing mean-zero assumption (k=0 mode) when inverting Δ / defining A^{-1}.
- Mis-stating the cancellation (needs divergence-free and projection).
- Wrong Sobolev exponents (3D embeddings; ∇u in L∞ needs s>5/2).
- Incorrect interpolation exponent (requires α>1).
- Failure to justify commutator estimate on T^3 (periodic setting).
- Confusion between A and -Δ (must include P and domain restrictions).
- Forcing regularity assumptions needed to propagate smoothness.

# Citation stubs (fill later)
- Lions: hyperdissipative NS global regularity for α ≥ (d+2)/4.
- Temam: NS functional analysis, local existence, energy methods.
- Kato–Ponce: commutator estimates.
- Tao: logarithmically supercritical regime (optional).
