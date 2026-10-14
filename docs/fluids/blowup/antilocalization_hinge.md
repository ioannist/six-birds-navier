# Anti-localization hinge for the Clay chain (NS-BU-30)

## 1) Wavepacket obstruction (no-go slot)

The mechanized result is a discrete diagonal-gain obstruction on a feasible
class containing coordinate basis tests. A continuum wavepacket-to-strain
return for the nonlinear NS bridge would additionally require an actual gain
representation, solenoidal feasible tests and a controlled fixed-background
limit. These inputs are OPEN; the existence of localized shell tests alone
does not establish that ICAP collapses to a BKM bound. See
`obstruction_certificate.md` for the exact scope.

## 2) Hinge statement (HL-P2-ANTILOC)

Use shell-scale Euclidean balls in `T^3`:

`B_j(x_0) := B(x_0, c 2^{-j})`.

Define the feasible shell-input set `U_j^feas`.

`HL-P2-ANTILOC:` for all `v in U_j^feas`,

`sup_{x_0 in T^3} int_{B_j(x_0)} |v(x)|^2 dx <= eta_j ||v||_2^2,  eta_j -> 0`.

## 3) Exclusion of the specific localized tests

If a proposed test sequence has at least a fixed positive fraction c0 of its
energy in one shell-scale ball, `HL-P2-ANTILOC` excludes those tests once
eta_j<c0. This blocks that particular proposed test argument. It does not
provide ICAP, the missing continuum no-go return, WORK, DIV or continuation.

## 4) Channelization implies anti-localization (addressability obstruction theorem)

Assume:

- `U_j^feas = Ran(G_j)` with `G_j : R^(m_j) -> U_j`,
- `v = sum_{r=1}^{m_j} a_r psi_{j,r}`,
- `{psi_{j,r}}_{r=1}^{m_j}` orthonormal in `L^2(T^3)`,
- channel delocalization bound:
  `sup_{x_0 in T^3} int_{B_j(x_0)} |psi_{j,r}(x)|^2 dx <= kappa 2^{-3j}` for each `r`.

Then for every `v in U_j^feas`,

`sup_{x_0 in T^3} int_{B_j(x_0)} |v(x)|^2 dx <= eta_j ||v||_2^2`,

with

`eta_j <= kappa m_j 2^{-3j}`.

### Proof

For `v = sum_r a_r psi_{j,r}`, Cauchy-Schwarz gives

`|v(x)|^2 <= (sum_r |a_r|^2)(sum_r |psi_{j,r}(x)|^2)`.

Integrate over `B_j(x_0)`:

`int_{B_j(x_0)} |v|^2 <= (sum_r |a_r|^2) sum_r int_{B_j(x_0)} |psi_{j,r}|^2`.

By orthonormality, `sum_r |a_r|^2 = ||v||_2^2`.
By delocalization of each channel, `int_{B_j(x_0)} |psi_{j,r}|^2 <= kappa 2^{-3j}`.
Therefore

`int_{B_j(x_0)} |v|^2 <= ||v||_2^2 (kappa m_j 2^{-3j})`.

Taking `sup_{x_0}` yields the claim.

For fixed kappa, the additional condition m_j 2^{-3j}→0 makes this returned
bound tend to zero. A finite channel count at each j alone does not supply
that asymptotic condition.

## 5) Measurement hooks (evidence only)

For a zero shell the local fraction is defined as zero in the diagnostics;
the analytic inequality itself needs no division by a zero norm. Numerical
windows are periodic grid cubes, whereas the analytic statement uses balls.
Comparing their constants requires the chosen size/geometry adapter.
Ranks are measured using uncentered snapshot squared energy, not centered
statistical variance. The discarded residual must be charged when returning
a top-m95 span bound to full snapshots; see section 6 of
channelized_feasibility_antilocalization.md.

These scripts provide evidence hooks for `eta_j`, `L_j`, `kappa_j`, and diagnostic `m95`, but they are not definitions or proofs:

- `python/scripts/nsbu28_antilocalization_diagnostic.py`
- `python/scripts/nsbu29_antilocalization_sweep.py`
