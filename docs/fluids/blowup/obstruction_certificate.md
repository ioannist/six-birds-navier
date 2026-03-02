# Setup + canonical objects
- LP shells S_j, blocks u_j := P_j u, low-pass u_{≤j}.
- Closure term C_j(u) and shell port output y_j := P_j C_j(u).
- Power p_j(t) := ⟨u_j(t), y_j(t)⟩_{L^2}, positive power p_j^+(t) := max(p_j(t),0).
- Positive work W_j^+[s,t] := ∫_s^t p_j^+(τ) dτ.

# What we need for No-Zeno (CAP/DIV target)
- Capacity form: W_j^+[s,t] ≤ Λ(j) ∫_s^t |u_j(τ)|_{L^2}^2 dτ.
- With budget ∫_s^t |u_j|_2^2 ≤ \bar B(j) (t-s), define Cap(j) := Λ(j)\bar B(j).
- No-Zeno engine (constant work quantum w(j)=θ):
  Δt_j ≥ θ/Cap(j) and ∑_j 1/Cap(j) = ∞ ⇒ t_j → ∞.

# Obstruction inequality and why it’s the bottleneck
- From LP/commutator structure (NS-BU-14):
  |p_j(t)| ≲ ||∇u(t)||_∞ · ||u_j(t)||_2^2.
- Any capacity bound derived from this route inherits the ||∇u||_∞ channel,
  or else pays an exponential Bernstein factor in j.

# Obstruction Certificate (informal theorem)
Any proof of a scale-uniform (or polynomial-in-j) positive-work capacity bound for the NS SGS port
that does not route through controlling a ||∇u||_∞-type quantity would constitute a new inequality
beyond the classical energy+Littlewood–Paley toolkit. Concretely: within standard LP calculus,
the best capacity bound reduces to the obstruction norm ||∇u||_∞ (or else yields an exponential 2^{β j}
prefactor), so certifying DIV from energy alone is impossible.

# Regularity link (classical bottleneck statement)
If one had a bound of the form ∫_0^T ||∇u(t)||_∞ dt < ∞, then standard Grönwall-type arguments
imply smoothness up to time T (hence no blow-up). Therefore any route that certifies No-Zeno by
bounding capacity must ultimately control a ||∇u||_∞ / vortex-stretching channel that is the
classical regularity bottleneck.

# Binary conclusion for our program
- To prove no blow-up via No-Zeno, we must prove a new scale-passivity / bounded dissipation-density
  lemma that yields slow growth of Λ(j) independent of ||∇u||_∞.
- Otherwise, this approach cannot settle Clay on its own (it restates the classical obstruction).

Discrete mechanization hook: `Fluids.icap_implies_pointwise_pos_bound_of_basis_tests`
(`formal/Fluids/IcapNoGo.lean`) formalizes the localized-test no-go pattern
(`basis tests in feasible set` + `ICAP` => pointwise positive gain bound).
