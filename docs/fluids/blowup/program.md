# Scope and non-claims
- Not a claim of solving Clay; this is a primitive-shaped reformulation target.
- Blow-up (if any) interpreted as closure breakdown / feasibility failure.

# Canonical multiscale decomposition
- Dyadic projections Δ_j or shell masks on |k| ~ 2^j.
- Shell energy: E_shell_j(t) := 1/2 ||u_j(t)||_2^2.
- Frontier defined via storage growth crossing times (see below).

# Core port/bridge objects (canonical pipeline)
- Port space U_j; signals u_j(t), y_j(t); bridge operator y_j = Z_j[u_j].
- Power p_j(t) := ⟨u_j(t), y_j(t)⟩_{U_j}.
- Positive supply p_j^+(t) := max(p_j(t), 0).
- Positive work W_j^+[s,t] := ∫_s^t p_j^+(τ) dτ.
- Storage S_j(x) ≥ 0 with passivity: S_j(t)-S_j(s) ≤ ∫ p_j ≤ W_j^+.

# Frontier and crossing times (storage-driven)
- Fix threshold θ > 0.
- A_{j+1}(t) := sup_{τ∈[t_j,t]} (S_j(τ) - S_j(t_j)).
- t_{j+1} := inf{t ≥ t_j : A_{j+1}(t) ≥ θ}.
- Frontier index: j_*(t) := max{j : t_j ≤ t}.

# Lemma chain (canonical)
1) passivity + storage-frontier ⇒ work quantum W_j^+[t_j,t_{j+1}] ≥ θ.
2) throughput + budget ⇒ W_j^+[t_j,t_{j+1}] ≤ Cap(j)·Δt_j.
3) hence Δt_j ≥ θ / Cap(j).
4) if ∑ 1/Cap(j) = ∞ ⇒ no Zeno (t_j → ∞).

# Hard open lemma frontier
- justify throughput bound Λ(j) (integrated, scale-local).
- justify budget bound \bar B(j).
- control growth of Cap(j) across scales.
- route mismatch control as a derivation tool.

# Measurement hooks (existing scripts)
- nsbu01_dyadic_frontier_demo.py
- nsbu02_latency_series.py
- nsbu03_bkm_proxy_vs_frontier.py
- nsbu04_sgs_frequency_response.py
- nsbu05_sgs_realization_fit.py
- nsbu06_route_dependence.py
