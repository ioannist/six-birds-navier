# NS packaged bridge object (Z_j)

## A) Setup and spaces (explicit)
- [ ] Domain: periodic T^3, divergence-free and mean-zero velocity fields.
- [ ] H := L^2_sigma(T^3) = {u in L^2(T^3; R^3) : div u = 0, mean(u)=0}.
- [ ] H^s_sigma(T^3) := H^s(T^3; R^3) ∩ H with s > 5/2.
- [ ] Fix a time horizon T > 0 and time-signal spaces:
  - L^2([0,T]; H), L^infty([0,T]; H), L^2([0,T]; H^{s+1}_sigma), L^infty([0,T]; H^s_sigma).
- [ ] Low-pass shell spaces:
  - U_{<=j} := {v in H : v = G_j v}, U_j := {v in H : v = Delta_j v}.
- [ ] Admissible port histories:
  - U^{adm}_{<=j}([0,T]) subset of L^infty([0,T]; H^s_sigma) ∩ L^2([0,T]; H^{s+1}_sigma) with G_j v = v and such that the high-pass IVP in Section E is well-posed and unique on [0,T].

## B) Dyadic decomposition and projectors
- [ ] Fix a Littlewood-Paley decomposition (Delta_j)_{j>=0} on T^3.
- [ ] Define G_j := sum_{ell<=j} Delta_ell, and P_{>j} := I - G_j.
- [ ] Decomposition: u_j := Delta_j u, u_{<=j} := G_j u, u_{>j} := P_{>j} u.
- [ ] Leray projector: P is the orthogonal projector in L^2 onto divergence-free fields.
- [ ] NS nonlinearity: NL(u) := P div(u ⊗ u).

## C) Dyadic equation for the shell and the environment influence term
- [ ] Full NS (projected): partial_t u = nu Delta u - NL(u).
- [ ] Apply Delta_j:
  - partial_t u_j = nu Delta u_j - Delta_j NL(u_{<=j}) + E_j(u_{<=j}, u_{>j}).
- [ ] Environment influence term (interactions involving u_{>j}):
  - E_j(u_{<=j}, u_{>j}) := -Delta_j NL(u_{<=j} + u_{>j}) + Delta_j NL(u_{<=j}).

## D) Closure / SGS term at level j
- [ ] Resolved closure on low modes:
  - C_{<=j}(u) := G_j NL(u) - G_j NL(u_{<=j}).
- [ ] Shell-level output used for capacity:
  - y_j(u) := Delta_j C_{<=j}(u).
- [ ] Port input is u_j := Delta_j u; port output is y_j(u) paired with u_j in power identities.

## E) Packaging map and definition of the bridge object (Z_j)
- [ ] Fix T > 0 and u_0 in H^s_sigma with s > 5/2 (strong solution regime).
- [ ] For a candidate resolved history v in U^{adm}_{<=j}([0,T]), define the environment state w as the unique solution (when it exists) of the high-pass equation:
  - partial_t w = nu Delta w - P_{>j} NL(v + w),  w(0) = P_{>j} u_0.
- [ ] Define the bridge output as the SGS/closure term induced by w:
  - (Z_j[v])(t) := y_j(v(t) + w(t)) = Delta_j C_{<=j}(v(t) + w(t)).
- [ ] Z_j is generally nonlinear, non-LTI, time-varying, and memoryful (depends on v through the IVP for w).
- [ ] Causality (explicit): if v^(1), v^(2) in U^{adm}_{<=j}([0,T]) agree a.e. on [0,t], then Z_j[v^(1)] and Z_j[v^(2)] agree a.e. on [0,t].
- [ ] Choice of object: Z_j is defined as a map on the admissible set U^{adm}_{<=j}([0,T]) where the environment IVP is well-posed and unique.

## F) One-paragraph boxed definition (required)

**Definition (NS packaged bridge object (Z_j) on ([0,T])).** Fix T>0, s>5/2, and u_0 in H^s_sigma(T^3). Let U^{adm}_{<=j}([0,T]) subset of L^infty([0,T]; H^s_sigma) ∩ L^2([0,T]; H^{s+1}_sigma) be the set of low-pass histories v with G_j v = v for which the high-pass IVP partial_t w = nu Delta w - P_{>j} NL(v+w), w(0)=P_{>j}u_0 admits a unique solution w in L^infty([0,T]; H^s_sigma) ∩ L^2([0,T]; H^{s+1}_sigma). The packaged bridge operator Z_j maps v in U^{adm}_{<=j}([0,T]) to Z_j[v] in L^2([0,T]; U_j) by Z_j[v](t) := Delta_j C_{<=j}(v(t)+w(t)), and is causal in the sense that agreement of two inputs a.e. on [0,t] implies agreement of outputs a.e. on [0,t].

## G) Diagnostics only (do not confuse with definition)
BU-04/05 (`python/scripts/nsbu04_sgs_frequency_response.py`, `python/scripts/nsbu05_sgs_realization_fit.py`) compute empirical frequency responses under strong simplifying assumptions (approximate LTI behavior, finite sampling, stationarity-ish windows). Those fitted H_s(omega) and positive-real realizations are diagnostics only and are not the definition of Z_j; they are optional measurement hooks that may approximate aspects of Z_j in regimes where LTI approximations appear reasonable.
