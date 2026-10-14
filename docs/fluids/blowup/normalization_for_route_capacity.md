# Normalization for route-capacity diagnostics

- U_rms^2 := time_avg(mean_x(|u^(0)(t,x)|^2)) over the same burn-in trimmed window.
- k_cut(n) := max{|k| : dealias_mask_3d(N_level(n)) is True} for each depth n.
- C_tilde(n) := C(n) / (U_rms * k_cut(n)).
- C_tilde is dimensionless and removes trivial amplitude scaling (u -> a u).
- U_rms * k_cut is an eddy-turnover rate built from energy, not from ||∇u||_∞.
- C_tilde is a diagnostic normalization only; it does not prove route-capacity hypotheses.
- Raw and normalized mismatch trends are reported separately. Decay after a
  depth-varying k_cut normalization does not establish raw summability.
- Finite trend PASS means only a sampled trend/tolerance check. Every finite
  prefix has both summable and nonsummable continuations. The exact telescoping
  identity, not a numerically small finite prefix, proves zero field mismatch
  for the nested projection model.
- A fixed 5% factor tolerance accumulates as 1.05^n, an exponential loss. It
  cannot be silently dropped from the polynomial-capacity theorem. Empty or
  all-zero capacity data have no fitted step ratio and are INCONCLUSIVE.

Interpretation:
For the spectral truncation ladder used in BU-21/23, route mismatch e_route is expected to be ~0
by telescoping/associativity, so the mismatch check is a sanity gate. The nontrivial hypothesis
for RouteCapacity is the step-factor condition on normalized capacities C_tilde.
