# Per atom ICAP and the full SGS source

Use the supplied-force convention in [ns_bridge_object.md](ns_bridge_object.md):
y_j=Δ_jG_j(NL(G_j u)-NL(u)). The exact shell energy balance is

    (1/2)d/dt||u_j||₂² + ν||∇u_j||₂²
      = -<u_j,Δ_j NL(G_j u)> + <u_j,y_j>.

## Conditional aggregation

If y_j=Σ_{r=1}^{m_j}y_{j,r} and every atom on the same actual feasible
history satisfies

    ∫_s^t max(<u_j,y_{j,r}>,0) ≤ Λ0 ∫_s^t ||u_j||₂²,

then max(Σ_r p_r,0)≤Σ_r max(p_r,0) gives total ICAP constant m_jΛ0.
The finite sum inequalities are mechanized in `Fluids.ECTCap`. The
continuum time integrals and existence of this NS atom decomposition are
separate, unmechanized inputs. Λ0 must be uniform in atoms, depth, histories
and the windows being claimed.

## What is not a source proof

The previous estimate with 2^{2j}||u_{>j}||₂||u_j||₂ as the complete
high–high remainder was not justified. The complete tensor is
v⊗w+w⊗v+w⊗w. Its verified finite-shell bound, including the quadratic
high-state contribution, is in [cap_bound_attempt.md](cap_bound_attempt.md).
The high–high triad there refutes the shell-only strain bound used by the
old attempt. No uniform per-atom ICAP supplier is proved here.

Storage passivity S(t)-S(s)≤∫p≤∫p⁺ does not bound p⁺ by shell energy.
One sufficient hypothesis is p_r⁺≤g_r||u_j||₂² with 0≤g_r≤Λ0 a.e.
Alternatively, ∫_s^t g_r≤Λ0(t-s) on **every** subinterval, with locally
integrable g_r, implies g_r≤Λ0 a.e. by Lebesgue differentiation. An
integral bound on just one window, or mere integrability, is insufficient
to multiply an arbitrary time-varying shell energy by Λ0. Feasible-history
quantification and the independent per-atom source remain OPEN.
