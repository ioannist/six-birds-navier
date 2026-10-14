# Target statement
- Not a Clay solution claim; this is a settlement map.
- Route A: No-Zeno latency lower bound ⇒ no finite-time infinite frontier.
- Route B: Blow-up ⇒ closure/passivity failure in SGS bridge.

# Setup and notation
- Domain: T^3; Leray projector P; Stokes operator A.
- Dyadic blocks Δ_j, low-pass S_j; u_j := Δ_j u.
- Vorticity ω := ∇×u (BKM linkage).
- Frontier via storage-growth crossing times (single definition).

# Ladder observables and “Zeno” definitions
- E_shell_j(t) := 1/2 ||u_j||_2^2.
- A_{j+1}(t) := sup_{τ∈[t_j,t]} (S_j(τ) - S_j(t_j)).
- t_j := inf{t : j_*(t) ≥ j}, Δt_j := t_{j+1}-t_j.
- Finite nonempty attained crossings and continuous storage/integrable power
  are required for the work-quantum proof; see scale_passivity_axioms.md.
- Zeno criterion: Σ Δt_j < ∞; No-Zeno: Σ Δt_j = ∞.
- Finite-resolution caveat: plateau ≠ theorem.

# Dyadic energy / flux identities (templates)
- Template: d/dt E_j + ν 2^{2j} E_j ≤ transfer + forcing.
- Flux schematic Π_j via dyadic sums (template only; not derived).

# Where incompressibility and vortex stretching enter
- Cancellation ⟨B(u,u),u⟩=0 uses ∇·u=0.
- Vorticity: ∂t ω + (u·∇)ω = (ω·∇)u + νΔω.
- (ω·∇)u is the hard BKM term.

# Candidate “scale-passivity” axioms (stronger than global energy)
- These are candidate premises/diagnostics, not established laws of the actual NS bridge.
- SP-1 Shellwise dissipativity: mean_t Σ Re(C·conj(u)) ≥ 0, with C in the
  extracted-power/LHS convention. The supplied closure has the opposite sign.
- SP-2 Sampled frequency-response diagnostic: Re(H_s(ω)) ≥ 0 for most sampled
  ω. This is not a positive-real/passivity theorem for the nonlinear bridge.
- SP-3 Throughput bound: high-k dissipation dominates beyond j0.
- SP-4 Route-robustness: hi→co ≈ hi→mid→co within ε(j).

# Lemma chain inventory (No-Zeno route)
- NZ-L0 Definitions: p, p^+, W^+, storage S, passivity inequality.
  - Claim: S_j(t)-S_j(s) ≤ ∫ p_j ≤ W_j^+.
- NZ-L1 Frontier definition via storage growth and crossing times.
  - Premise: every required level has a nonempty hitting set and a finite attained
    crossing. Existence cannot be concluded just from writing an infimum.
- NZ-L2 Passivity ⇒ work quantum.
  - Conditional continuous-time claim: W_j^+[t_j,t_{j+1}] ≥ θ, using passivity and
    attained threshold. This is not the theorem proved by the division-only Lean file.
- NZ-L3 Throughput + budget assumptions.
  - Claim: W_j^+[t_j,t_{j+1}] ≤ Cap(j)·Δt_j.
- NZ-L4 Latency lower bound.
  - Claim: Δt_j ≥ θ / Cap(j).
- NZ-L5 Divergence ⇒ No-Zeno.
  - Claim: if ∑ 1/Cap(j)=∞ then t_j → ∞ (Lean: Fluids.Zeno).
- NZ-L6 Route mismatch ⇒ capacity growth control (slot).
  - Claim: Cap growth bounded by mismatch functional.

# Lemma chain inventory (Blow-up ⇒ closure/passivity failure route)
- This is an unproved proposed route; BF-L1 through BF-L4 are not consequences
  of the abstract No-Zeno results or of the frequency-response diagnostics.
- BF-L0 Hypothetical finite terminal time with unbounded continuation norm.
- BF-L1 OPEN: the actual chosen storage frontier must cross arbitrarily high
  levels before T. LP smoothness at each fixed time does not supply this return.
- BF-L2 OPEN: a quantified actual-path passivity failure would need its own
  argument; blow-up alone has not been shown to force it here.
- BF-L3 OPEN: a nontrivial packaged-history route mismatch would need its own
  estimate; the nested instantaneous field-level composition has zero mismatch.
- BF-L4 OPEN: a mathematical failure criterion for P6/P4 feasibility is required.
- BF-L5 Diagnostics: BU-04/05/06 detect failure modes.
- BF-L6 Identify minimal SP-* that fails first.

# Proof stress-test checklist (exponent / embedding / bookkeeping pitfalls)
- Sign mistakes: using p instead of p^+.
- Forgetting storage baseline at t_j.
- Mixing bridge mismatch with restriction idempotence.
- A fast-growing upper capacity bound may fail to certify divergence; it does
  not prove convergence of the actual reciprocal series.
- Wrong scaling exponents in 3D (criticality).
- Confusion between energy/enstrophy/Sobolev norms.
- Mean-zero requirement for A^{-1}.
- Dyadic vs sharp spectral shells differences.
- Bony decomposition index bookkeeping.
- De-aliasing / finite-grid artifacts (numerics).
- Forcing invalidates naive energy arguments.
- Shellwise vs global passivity mismatch (backscatter).

# Citation stubs (fill later)
- [REF: Leray]
- [REF: CKN]
- [REF: BKM]
- [REF: Kato–Ponce]
- [REF: Littlewood–Paley]
- [REF: Lions hyperviscous]
