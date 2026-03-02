# Channelized feasibility and anti-localization (NS-BU-28)

## Scope and intent

This note defines a feasibility-restricted input class for shell-level tests and records the exact anti-localization lemma slot needed by the NS blow-up program.

This is a math-level hinge statement, not a Clay proof.

## 1) Channelized feasibility class (port language)

Fix dyadic shell index `j` and define

`U_j := Ran(Delta_j) subset L^2(T^3)` (optionally intersected with divergence-free fields).

Define feasible shell inputs as the range of a synthesis map

`G_j : R^(m_j) -> U_j`,

`u = G_j a = sum_{r=1}^{m_j} a_r psi_{j,r}`,

where `psi_{j,r} in U_j` are channel shapes and `m_j` is the channel count.

Feasible-input ICAP statements are then understood as quantifying only over `u in Ran(G_j)`, not all of `U_j`.

## 2) Delocalization condition at shell scale

Let `B(x_0, c 2^{-j})` denote a Euclidean ball in `T^3` with radius `c 2^{-j}`.

Assume each channel shape satisfies the scale-local mass bound

`sup_{x_0 in T^3} int_{B(x_0, c 2^{-j})} |psi_{j,r}(x)|^2 dx <= kappa 2^{-3j}`.

Interpretation: each addressable channel is spread across many uncertainty-scale cells, with no single shell-scale ball carrying large `L^2` mass.

## 3) Anti-localization lemma (orthonormal channels)

Assume `{psi_{j,r}}_{r=1}^{m_j}` is orthonormal in `L^2(T^3)` and obeys the delocalization condition above. Then for every `u in Ran(G_j)`,

`sup_{x_0 in T^3} int_{B(x_0, c 2^{-j})} |u(x)|^2 dx <= eta_j ||u||_2^2`,

with

`eta_j := kappa m_j 2^{-3j}`.

### Proof

Write `u = sum_{r=1}^{m_j} a_r psi_{j,r}`. Pointwise,

`|u(x)|^2 = |sum_r a_r psi_{j,r}(x)|^2 <= (sum_r |a_r|^2)(sum_r |psi_{j,r}(x)|^2)`

by Cauchy-Schwarz in `R^(m_j)` (or `C^(m_j)`).

Integrating over any shell-scale ball `B = B(x_0, c 2^{-j})` gives

`int_B |u|^2 <= (sum_r |a_r|^2) sum_r int_B |psi_{j,r}|^2`.

By orthonormality, `sum_r |a_r|^2 = ||u||_2^2`.
By the channel delocalization hypothesis, `int_B |psi_{j,r}|^2 <= kappa 2^{-3j}` for each `r`.
Hence

`int_B |u|^2 <= ||u||_2^2 * m_j * kappa * 2^{-3j}`.

Taking supremum over `x_0` yields the claim.

## 4) Explicit link to the wavepacket obstruction

If `eta_j -> 0` as `j -> infty`, then feasible inputs in `Ran(G_j)` cannot realize shell-scale wavepacket localization.
Therefore, ICAP posed only on feasible inputs does not automatically imply the classical `L^infty` strain control route based on wavepacket test functions, because those test functions are not feasible-addressable in the packaged interface language.

## 5) Lemma slot for the Clay chain

Use this as the placeholder lemma:

- `HL-P2-ANTILOC`: packaged channelized feasibility implies anti-localization factor `eta_j`.
- Bridge use: replace "ICAP for all shell tests" by "ICAP on feasible tests + anti-localization."
- Decision hinge: whether NS dynamics admits such channelized feasible classes with `eta_j -> 0` at depth.

Lean mechanization: the finite discrete version is formalized as
`Fluids.mass_cell_le_card_mul_deloc_mul_coeffEnergy` in `formal/Fluids/AntiLocalization.lean`.
