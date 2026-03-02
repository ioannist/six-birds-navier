# NS-BU-25 Clay theorem chain scaffold (paper-grade)

## A. What this repository does and does not claim
- This repository mechanizes abstract No-Zeno/route-capacity algebra and records NS-facing hinge lemmas and diagnostics.
- This repository does **not** claim a completed proof of the Clay 3D Navier–Stokes regularity theorem.
- Any Clay-facing conclusion in this document is conditional on explicit OPEN assumptions marked below.

## B. Two branches of the chain (explicit split)
- **Branch 1 (no-go for classical Clay via ICAP certificates):**
  ICAP over localization-rich feasible tests implies a pointwise positive-gain bound in the discrete model (`formal/Fluids/IcapNoGo.lean`), mirroring the wavepacket obstruction in `docs/fluids/blowup/obstruction_certificate.md`.
- **Branch 2 (SBT-legal NS variant, conditional):**
  If feasible inputs are structurally anti-localized (HL-P2-ANTILOC / HL-P5-CHANNEL-DEL, OPEN) and persistence assumptions hold, the mechanized No-Zeno engine applies to that restricted class.
- Every hinge not mechanized for NS itself is tagged **ASSUMPTION / OPEN**; these are settlement points, not proved facts here.

## C. Target PDE statement (classical 3D NS, Clay-compatible)
- [ ] Domain fixed: 𝕋^3 (periodic). Remark: the same chain can be stated on ℝ^3 with the usual decay assumptions.
- [ ] PDE (ν > 0):
  - ∂_t u + (u·∇)u + ∇p = ν Δu,  ∇·u = 0,  u(0) = u_0.
- [ ] Data: u_0 is divergence-free and smooth (or u_0 ∈ H^k with k > 5/2).
- [ ] Target claim (paper theorem form):
  - For every divergence-free u_0 ∈ H^k(𝕋^3), k > 5/2, the unique strong solution is global and smooth (no finite-time blow-up).

## D. Canonical multiscale decomposition
- [ ] Dyadic shells: S_j := {k ∈ ℤ^3 : 2^j ≤ |k| < 2^{j+1}}.
- [ ] Projectors: P_j = Fourier projector to S_j, G_j = low-pass to |k| < 2^{j+1}.
- [ ] Fields: u_j := P_j u, u_{≤j} := G_j u, P_{≤j} := G_j.
- [ ] Shell energy: E_shell_j(t) := (1/2) ||u_j(t)||_2^2.

## E. Canonical SGS “port” objects in NS terms
- [ ] SGS closure at cutoff j:
  - C_j(u) := P_{≤j} 𝒫 ∇·(u_{≤j}⊗u_{≤j} − (u⊗u)_{≤j}), with 𝒫 the Leray projector.
- [ ] Port output (shell-localized): y_j(t) := P_j C_j(u(t)).
- [ ] Port input: u_j(t) := P_j u(t).
- [ ] Power and positive work:
  - p_j(t) := ⟨u_j(t), y_j(t)⟩_{L^2}.
  - p_j^+(t) := max(p_j(t), 0).
  - W_j^+[s,t] := ∫_s^t p_j^+(τ) dτ.

## F. Capacity definition (ICAP / positive-work capacity)
- [ ] ICAP (positive-work form): for all [s,t] ⊂ [0,T],
  - W_j^+[s,t] ≤ Λ(j) ∫_s^t ||u_j(τ)||_2^2 dτ.
- [ ] Optional throughput budget (if using Cap(j)):
  - ∫_s^t ||u_j(τ)||_2^2 dτ ≤ B(j) (t-s) for all [s,t].
  - Cap(j) := Λ(j) B(j).
- [ ] Interpretation: Λ(j) is the time-integrated dissipation density / bandwidth of the SGS port.

## G. Storage frontier (work-quantized crossing rule)
- [ ] Storage: S_j(t) ≥ 0 for eliminated scales u_{>j}.
- [ ] Activity:
  - A_{j+1}(t) := sup_{τ∈[t_j,t]} (S_j(τ) − S_j(t_j)).
- [ ] Crossing times:
  - t_{j+1} := inf{ t ≥ t_j : A_{j+1}(t) ≥ θ }.
- [ ] Work quantum: w(j) := θ (constant, canonical).
- [ ] WORK quantum lemma (slot):
  - W_j^+[t_j, t_{j+1}] ≥ θ.

## H. ECT block (capacity growth control)
- [ ] **Theorem (ECT ⇒ linear capacity growth).**
  - Hypotheses: the SGS port at depth j decomposes into m_j atoms (parallel sum) with
    - (i) m_j ≤ C0 (j+1) (mode/sector compression),
    - (ii) each atom satisfies ICAP with the same constant Λ0.
  - Conclusion: the aggregated SGS port satisfies ICAP with
    - Λ(j) ≤ Λ0 * m_j ≤ Λ0 * C0 * (j+1).
- [ ] One-line inequality used:
  - max(∑_r a_r, 0) ≤ ∑_r max(a_r, 0) (pointwise in time).
- [ ] **HL-P2-ANTILOC (new hinge between ECT/ICAP and wavepacket tests).**
  - Feasible shell inputs satisfy:
    - sup_{x0} ∫_{B(x0,c2^{-j})} |v|^2 ≤ η_j ||v||_2^2, with η_j → 0.
  - Use in chain:
    - Any step that would quantify over wavepacket-like feasible shell tests must be replaced by feasibility-restricted anti-localized tests.
  - Why this is the new hinge:
    - Uniform ICAP over localization-rich sets collapses to a BKM-type channel; anti-localization of feasible inputs is the only coherent bypass currently identified.

## I. No-Zeno theorem chain (Lean references)
1) [ ] **WORK + CAP ⇒ latency bound.**
   - Δt_j := t_{j+1} − t_j satisfies Δt_j ≥ w(j)/Cap(j) (or Δt_j ≥ θ/(B(j)Λ(j))).
   - Reference: `formal/Fluids/ZenoWorkCap.lean`.
2) [ ] **Divergence ⇒ no finite upper bound on crossing times.**
   - If ∑_j w(j)/Cap(j) = ∞ then no finite T with t_j ≤ T for all j.
   - Reference: `formal/Fluids/Zeno.lean`, `formal/Fluids/ZenoWorkCap.lean`.
3) [ ] **Linear Λ(j) ⇒ divergence (analytic lemma slot).**
   - If Cap(j) ≤ K (j+1) for large j, then ∑_j 1/Cap(j) = ∞.
   - Status: standard analytic fact (harmonic divergence) or future Lean mechanization ticket.

## J. Clay bridge lemma slot (No-Zeno ⇒ no blow-up)
- [ ] Define dyadic H^k tail for k > 5/2:
  - Tail_k(J,t) := ∑_{ℓ≥J} 2^{2kℓ} ||u_ℓ(t)||_2^2.
- [ ] **HL-NS-ZENO (Blow-up ⇒ Zeno frontier).**
  - If a strong solution blows up at time T* in H^k, then t_j ↑ T* (infinitely many frontier crossings in finite time).
  - Mechanism needed: sup_{t<T*} ||u(t)||_{H^k} = ∞ ⇒ Tail_k(J,t) → ∞ as J→∞ ⇒ forces infinitely many storage frontier crossings.
- [ ] This is the single PDE-hard bridge lemma that turns No-Zeno into Clay regularity.

## K. Final conditional theorem (paper-grade)
- [ ] **Theorem (Conditional Clay regularity from ECT + Zeno bridge).**
  - Assume:
    1) NS admits the packaged SGS port decomposition satisfying the ECT hypotheses uniformly in j,
    2) feasible shell inputs satisfy HL-P2-ANTILOC (anti-localization, **ASSUMPTION / OPEN**),
    3) the storage frontier satisfies WORK-quantization with w(j)=θ,
    4) throughput holds for the chosen port input (defines Cap(j)), and
    5) HL-NS-ZENO holds (**ASSUMPTION / OPEN**).
  - Then there is no finite-time blow-up of the strong solution; hence the 3D incompressible NS solution is global and smooth.

## L. Assumption audit and binary settlement points
- [ ] Mechanized (Lean):
  - Zeno logic and divergence criterion: `formal/Fluids/Zeno.lean`.
  - Work/Cap ⇒ latency bound: `formal/Fluids/ZenoWorkCap.lean`.
  - Route capacity lemma: `formal/Fluids/RouteCapacity.lean`.
  - ECT algebraic step (finite sums): `formal/Fluids/ECTCap.lean`.
- [ ] Validated numerically (Python):
  - Capacity estimation: `python/src/capacity_estimator.py`.
  - Route diagnostics: `python/scripts/nsbu21_route_capacity_numbers.py`, `python/scripts/nsbu23_route_capacity_hypothesis_check.py`.
- [ ] Open and Clay-hard:
  - NS instantiation of ECT hypotheses (mode compression + per-atom ICAP uniformity).
  - NS proof of HL-P2-ANTILOC (or sufficient channelized feasibility theorem implying η_j → 0).
  - HL-NS-ZENO bridge (blow-up ⇒ Zeno frontier).
- [ ] Binary settlement points:
  - If ECT-for-NS and HL-NS-ZENO hold, the Clay regularity claim follows via the No-Zeno chain.
  - If either fails, the chain breaks exactly at that lemma slot (settlement frontier is explicit).

## M. Measurement hooks (repo links)
- [ ] Dyadic definitions: `python/src/nswave/dyadic.py`.
- [ ] Latency series: `python/scripts/nsbu02_latency_series.py`.
- [ ] BKM proxy overlay: `python/scripts/nsbu03_bkm_proxy_vs_frontier.py`.
- [ ] SGS frequency response: `python/scripts/nsbu04_sgs_frequency_response.py`.
- [ ] PR fit: `python/scripts/nsbu05_sgs_realization_fit.py`.
- [ ] Route capacity: `python/scripts/nsbu21_route_capacity_numbers.py`.
- [ ] Hypothesis check: `python/scripts/nsbu23_route_capacity_hypothesis_check.py`.
