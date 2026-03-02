# Normalization for route-capacity diagnostics

- U_rms^2 := time_avg(mean_x(|u^(0)(t,x)|^2)) over the same burn-in trimmed window.
- k_cut(n) := max{|k| : dealias_mask_3d(N_level(n)) is True} for each depth n.
- C_tilde(n) := C(n) / (U_rms * k_cut(n)).
- C_tilde is dimensionless and removes trivial amplitude scaling (u -> a u).
- U_rms * k_cut is an eddy-turnover rate built from energy, not from ||∇u||_∞.
- C_tilde is a diagnostic normalization only; it does not prove route-capacity hypotheses.

Interpretation:
For the spectral truncation ladder used in BU-21/23, route mismatch e_route is expected to be ~0
by telescoping/associativity, so the mismatch check is a sanity gate. The nontrivial hypothesis
for RouteCapacity is the step-factor condition on normalized capacities C_tilde.
