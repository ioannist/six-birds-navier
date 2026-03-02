# ECT H3 attempt: per-atom ICAP from NS structure

## A) Setup (exact objects)
- Domain: Omega = T^3, s > 5/2, time horizon [0,T].
- Dyadic blocks: u_j := Delta_j u, u_{<=j} := G_j u, P_{>j} := I - G_j.
- Leray projector: P onto divergence-free fields.
- NS nonlinearity: NL(u) := P div(u ⊗ u).
- Closure / SGS term (low-pass):
  - C_{<=j}(u) := G_j NL(u) - G_j NL(u_{<=j}).
  - Shell output: y_j(u) := Delta_j C_{<=j}(u).
- Power and positive work:
  - p_j(t) := <u_j(t), y_j(t)>_{L^2}.
  - p_j^+(t) := max(p_j(t), 0).
  - W_j^+[s,t] := int_s^t p_j^+(tau) d tau.

## B) What “atoms” mean in this context
- Hypothetical atomic/sector decomposition at level j:
  - y_j(t) = sum_{r=1}^{m_j} y_{j,r}(t).
  - p_{j,r}(t) := <u_j(t), y_{j,r}(t)>_{L^2}.
  - W_{j,r}^+[s,t] := int_s^t max(p_{j,r}(tau), 0) d tau.
- Target ICAP-per-atom inequality (H3):
  - W_{j,r}^+[s,t] <= Lambda_0 * int_s^t ||u_j(tau)||_2^2 d tau,
    with Lambda_0 independent of j and r.

## C) Baseline bound from dyadic energy identity (no new SBT hinge)
- Shell energy identity (exact):
  - (1/2) d/dt ||u_j||_2^2 + nu ||nabla u_j||_2^2 = -<u_j, Delta_j NL(u)>.
- Use y_j = Delta_j(NL(u) - NL(u_{<=j})) and write
  - <u_j, y_j> = <u_j, Delta_j NL(u)> - <u_j, Delta_j NL(u_{<=j})>.
- Expand NL(u) - NL(u_{<=j}) via u = u_{<=j} + u_{>j}:
  - u ⊗ u - u_{<=j} ⊗ u_{<=j} = u_{<=j} ⊗ u_{>j} + u_{>j} ⊗ u_{<=j} + u_{>j} ⊗ u_{>j}.
- Standard commutator/LP estimates (toolkit-only) give the obstruction channel:
  - |<u_j, y_j>| <= C ||nabla u||_infty * ||u_j||_2^2
    + C 2^{2j} ||u_{>j}||_2 * ||u_j||_2
    (dimensionless C, no dependence on j beyond the explicit dyadic factor).
- Energy-only fallback (using ||u_{>j}||_2 <= ||u||_2):
  - |<u_j, y_j>| <= C ||nabla u||_infty * ||u_j||_2^2 + C 2^{2j} ||u||_2 * ||u_j||_2.
- Conclusion of the baseline toolkit bound:
  - Any ICAP bound from this route must control ||nabla u||_infty (or pay exponential-in-j growth).

## D) Attempt the ECT route for per-atom ICAP
- Windowed (causal) ICAP semantics:
  - For any [s,t] subset [0,T], ICAP is stated as a windowed inequality with no reliance on earlier history,
    W_{j,r}^+[s,t] <= Lambda_0 * int_s^t ||u_j(tau)||_2^2 d tau.
- Sufficient condition template (integrated, per-atom):
  - p_{j,r}^+(tau) <= g_{j,r}(tau) * ||u_j(tau)||_2^2 a.e., and
    int_s^t g_{j,r}(tau) d tau <= Lambda_0 (t-s) for all [s,t].
- Why passivity is insufficient:
  - Passivity gives p_{j,r} <= W_{j,r}^+ but does not bound p_{j,r}^+ by a uniform multiple of ||u_j||_2^2.
  - A bounded-dissipation-density property is needed (lemma slot):
    "For each atom r, the positive power density is uniformly dominated by a time-integrable scalar bound
    independent of (j,r)."
- In NS, the only toolkit scalar control available for the SGS coupling is ||nabla u||_infty or a dyadic
  growth factor from high-high interactions; neither is uniform in j without extra structure.

## E) Final obstruction line (failure case)
W_{j,r}^+[s,t] <= int_s^t G(tau) * ||u_j(tau)||_2^2 d tau,  with G(tau) = ||nabla u(tau)||_infty (or, absent this control, an explicit dyadic-growth proxy from high-high interactions).
Uniform ICAP would require controlling G independently of (j,r), which is exactly the classical bottleneck channel.
