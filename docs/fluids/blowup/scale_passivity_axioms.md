# Scope and non-claims
- This document defines a math-only “decision engine” for whether a multiscale system can reach infinite ladder depth in finite time (Zeno cascade).
- This is not a PDE proof; it isolates the exact lemma frontier.

# Core objects
- Port space U_j (real Hilbert space).
- Port signals u_j(t) in U_j (input) and y_j(t) in U_j (response).
- Bridge operator Z_j: causal operator with y_j = Z_j[u_j] (not assumed linear).
- This single-port model is abstract. The actual NS bridge in
  `ns_bridge_object.md` takes the entire resolved history v and pairs its
  output with Δ_j v; it is not proved to factor through the shell history alone.
- Instantaneous power: p_j(t) := ⟨u_j(t), y_j(t)⟩_{U_j}.
- Positive supply: p_j^+(t) := max(p_j(t), 0).
- Positive supplied work: W_j^+[s,t] := ∫_s^t p_j^+(τ) dτ.
- Note: p^+ avoids cancellation between positive and negative power.

# Storage and passivity axiom
- Abstract state x_j(t) and storage S_j(x) ≥ 0.
- Scale-passivity / accounting axiom:
  S_j(x_j(t)) - S_j(x_j(s)) ≤ ∫_s^t p_j(τ) dτ ≤ W_j^+[s,t].

# Frontier and crossing times
- Fix threshold θ > 0.
- Define storage-growth activity:
  A_{j+1}(t) := sup_{τ∈[t_j,t]} (S_j(x_j(τ)) - S_j(x_j(t_j))) for t ≥ t_j.
- Crossing time:
  t_{j+1} := inf{t ≥ t_j : A_{j+1}(t) ≥ θ}.
- For this construction assume S_j is continuous, p_j is locally integrable,
  and the hitting set is nonempty with a finite crossing inside the time domain.
  Then the running maximum is continuous and the crossing threshold is attained.
- Frontier index: j_*(t) := sup{j : t_j ≤ t} in ℕ∪{∞}; a finite maximum is
  available only when finitely many crossings have occurred. It cannot be
  assumed in an argument whose purpose is to rule out Zeno accumulation.

# Lemma: frontier advance forces a work quantum
- If t_{j+1} is defined and passivity holds, then W_j^+[t_j,t_{j+1}] ≥ θ.
- Proof: continuity gives A_{j+1}(t_{j+1})≥θ. Compactness of [t_j,t_{j+1}]
  gives τ* attaining that maximum. Then θ≤S_j(τ*)-S_j(t_j)
  ≤W_j^+[t_j,τ*]≤W_j^+[t_j,t_{j+1}]. An infimum without endpoint
  regularity or an independently supplied threshold-attainment premise is not
  this proof. This continuous storage argument is not formalized by ZenoWorkCap.

# Throughput and feasibility (capacity) assumptions
- Integrated throughput bound: ∃Λ(j) > 0 s.t. for any [s,t] ⊂ [t_j,t_{j+1}],
  W_j^+[s,t] ≤ Λ(j) ∫_s^t |u_j(τ)|^2 dτ.
- Port feasibility budget: ∃\bar B(j) > 0 s.t. |u_j(t)|^2 ≤ \bar B(j) for a.e. t in [t_j,t_{j+1}].
- Capacity per unit time: Cap(j) := Λ(j) \bar B(j).

# Lemma: capacity ⇒ latency lower bound
- From throughput + budget:
  W_j^+[t_j,t_{j+1}] ≤ Cap(j) (t_{j+1}-t_j).
- Combine with work quantum to get:
  Δt_j := t_{j+1}-t_j ≥ θ / Cap(j).

# No-Zeno criterion via divergence
- If ∑_{j≥j0} θ / Cap(j) = ∞ (equivalently ∑ 1/Cap(j) diverges), then t_j → ∞ and no Zeno cascade.
- Sufficient conditions (checkable):
  - If Cap(j) ≤ C (j+1)^α for large j with α ≤ 1, then ∑ 1/Cap(j) diverges.
  - If Cap(j) has an exponential lower bound, the reciprocal sum converges
    (Zeno not ruled out by this certificate). An upper bound alone does not decide DIV.

# Route mismatch lemma slot (for later)
- Packaging operator Pack_{j←k} yields bridge at level j after eliminating >j.
- Mismatch scalar: RM(j) := ||Pack_{j←j+2} - (Pack_{j←j+1} ∘ Pack_{j+1←j+2})||.
- In a separate bounded-linear-map model one can set Cap(Z):=||Z||_{L²→L²}.
  The nonlinear NS bridge has no linear operator norm by default; its positive-work
  capacity is a different functional. A composition in this template requires
  matching domains/codomains and a proof of the bound in the currency actually used.
- Lemma slot (not proven here):
  Cap(Pack_{j←j+2}) ≤ Cap(Pack_{j←j+1}) · Cap(Pack_{j+1←j+2}) + RM(j).

# Measurement hooks (tie to existing repo artifacts)
- frontier/heatmap: python/scripts/nsbu01_dyadic_frontier_demo.py
- latency extraction: python/scripts/nsbu02_latency_series.py
- classical proxy overlay: python/scripts/nsbu03_bkm_proxy_vs_frontier.py
- SGS frequency response: python/scripts/nsbu04_sgs_frequency_response.py
- passive rational fit attempt: python/scripts/nsbu05_sgs_realization_fit.py
- route dependence diagnostic: python/scripts/nsbu06_route_dependence.py

# Mechanization link
- formal/Fluids/Zeno.lean contains the abstract “increments ≥ unbounded partial sums ⇒ no finite upper bound” skeleton.
