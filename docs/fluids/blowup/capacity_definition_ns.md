# Scope and non-claims
- This document defines a math-only decision engine for scale-frontier reachability in finite time.
- This is not a PDE proof; passivity/throughput are treated as axioms/templates until proven.

# State, decomposition, and shells
- Domain: T^3, u(x,t) ∈ ℝ^3, ∇·u=0.
- Fourier: û(k,t) with Parseval; dyadic shells S_j = {k: 2^j ≤ |k| < 2^(j+1)}.
- Shell projection: P_j u (Fourier mask), shell energy via `python/src/nswave/dyadic.py`.

# Canonical coarse-graining map
- Projection G_j keeps |k| < 2^(j+1); u_{≤j} := G_j u, u_{>j} := u - u_{≤j}.
- Write P_{≤j} := G_j for the low-pass projector.
- Nested projections: G_j G_{j+1} = G_j.

# Port variables (input/output)
- Port input: u_j(t) := P_j u(t) in U_j := L^2(T^3;ℝ^3) restricted to shell j.
- Leray projector: P_Leray onto divergence-free vector fields on T^3.
- SGS closure term:
  C_j(u) := P_{≤j} P_Leray ∇·(u_{≤j}⊗u_{≤j} - (u⊗u)_{≤j}).
- Port output: y_j(t) := P_j C_j(u(t)).
- Sign convention: C_j appears on RHS (closure forcing).

# Positive supplied work (W_j^+)
- Power: p_j(t) := ⟨u_j(t), y_j(t)⟩_{L^2}.
- Positive power: p_j^+(t) := max(p_j(t), 0).
- Positive work: W_j^+[s,t] := ∫_s^t p_j^+(τ) dτ.
- Use W^+ to avoid cancellation between positive/negative power.

# Capacity (Λ(j)) (integrated throughput bound)
- Capacity constant: smallest Λ(j) such that
  ∫_s^t ⟨u_j, y_j⟩ dτ ≤ Λ(j) ∫_s^t |u_j(τ)|_{L^2}^2 dτ for all [s,t].
- Positive-work version used for latency:
  W_j^+[s,t] ≤ Λ(j) ∫_s^t |u_j(τ)|_{L^2}^2 dτ. (This Λ(j) is the positive-work capacity used in the Zeno engine.)
- Budget (time-average): ∫_s^t |u_j(τ)|_{L^2}^2 dτ ≤ \bar B(j) (t-s).
- Per-crossing capacity: Cap(j) := Λ(j) \bar B(j), so W_j^+[s,t] ≤ Cap(j) (t-s).

# Work quantum (w(j)) from storage-frontier
- Abstract storage S_j(t) ≥ 0 for eliminated scales u_{>j}.
- Passivity template:
  S_j(t) - S_j(s) ≤ ∫_s^t ⟨u_j, y_j⟩ dτ ≤ W_j^+[s,t].
- Storage activity: A_{j+1}(t) := sup_{τ∈[t_j,t]} (S_j(τ) - S_j(t_j)).
- Crossing time: t_{j+1} := inf{t ≥ t_j : A_{j+1}(t) ≥ θ}.
- Work quantum: w(j) := θ, and W_j^+[t_j,t_{j+1}] ≥ θ.

# Non-rescalability / gauge invariance check
- Rescale: u_j' = a_j u_j, y_j' = a_j^{-1} y_j ⇒ p_j' = p_j and W_j^+ invariant.
- Capacity inequality changes under rescaling unless the port norm is fixed.
- Canonical choice: L^2 norm on u_j; Cap(j) is defined relative to this fixed Hilbert structure.

# How this connects to existing repo scripts
- nsbu01_dyadic_frontier_demo.py: shell energies, frontier proxies.
- nsbu03_bkm_proxy_vs_frontier.py: ω_∞ proxy + frontier overlays.
- nsbu04_sgs_frequency_response.py: estimates H_s(ω) from τ and u_{≤j}.
- nsbu05_sgs_realization_fit.py: fits positive-real rational families.
- nsbu06_route_dependence.py: route mismatch diagnostic.
- This document standardizes the math definitions these scripts approximate.

# Risk notes / failure modes
- Capacity Λ(j) may depend on ||∇u||_∞ (classical obstruction).
- Storage S_j with passivity inequality may not exist for NS.
- Port choice u_j,y_j may be too weak to control inter-scale transfer.
