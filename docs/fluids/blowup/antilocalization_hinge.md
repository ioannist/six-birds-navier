# Anti-localization hinge for the Clay chain (NS-BU-30)

## 1) Wavepacket obstruction (no-go slot)

If an ICAP-type inequality is required on a feasible class that still contains shell-scale wavepacket-localizing tests in shell `j`, then the usual wavepacket test route can recover an `L^infty`-strain-type bottleneck (BKM channel). In that setting, ICAP is not a bypass; it collapses back to the classical obstruction mechanism.

## 2) Hinge statement (HL-P2-ANTILOC)

Use shell-scale Euclidean balls in `T^3`:

`B_j(x_0) := B(x_0, c 2^{-j})`.

Define the feasible shell-input set `U_j^feas`.

`HL-P2-ANTILOC:` for all `v in U_j^feas`,

`sup_{x_0 in T^3} int_{B_j(x_0)} |v(x)|^2 dx <= eta_j ||v||_2^2,  eta_j -> 0`.

## 3) Why this kills the wavepacket obstruction

The no-go route needs localization-rich feasible tests (wavepacket-like shell inputs) to convert feasible ICAP into local strain control. `HL-P2-ANTILOC` directly denies that richness: shell-scale local mass is uniformly suppressed by `eta_j -> 0`, so the wavepacket test class is not feasible-addressable.

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

## 5) Measurement hooks (evidence only)

These scripts provide evidence hooks for `eta_j`, `L_j`, `kappa_j`, and diagnostic `m95`, but they are not definitions or proofs:

- `python/scripts/nsbu28_antilocalization_diagnostic.py`
- `python/scripts/nsbu29_antilocalization_sweep.py`
