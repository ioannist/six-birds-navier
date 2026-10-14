# Completion audit of the original Navier Stokes scaffold

The requested scope is the mathematics and implementation supporting
*Six Birds for Navier–Stokes: A Mechanized No-Zeno Scaffold*. This audit
covers that manuscript's theorem chain, the nine original Lean modules,
the associated analytic dossiers and Python experiment implementation.
It does not resume or claim to review the separate unconditional campaign,
the later sibling Navier paper, or historical archived implementations.

## Requirement and evidence

| Requirement | Result and evidence |
| --- | --- |
| Preserve the starting state | Commit `b2c2316` precedes the mathematical repairs. Original reports and rewritten-dossier snapshots are preserved. |
| Thorough mathematical review | [Review](mathematics_review.md) maps each original abstract result to its actual premises and conclusion, identifies the false SGS estimate by an actual solenoidal datum, and corrects projection, sign, crossing and reciprocal-comparison arguments. |
| Preserve claim strength when justified | Main abstract theorem statements and title retain their strength. Fixed-μ hyperviscous regularity at α≥5/4 is repaired, and the existing conditional smooth-interval recovery is restored with an explicit O(μ) rate. Scalar-symbol recovery is conditional on actual neighborhood passivity and time locality. |
| Repair mechanizations | Actual tame-factor, orthonormal Parseval, field-energy, feasible-basis and derivative-operator adapters are proved. The fixed-mode counterexample is mechanized. Old exports remain available. |
| Check actual Lean dependencies | `lake build` and `lake env lean ScaffoldAxiomAudit.lean` pass with Lean/mathlib 4.8.0. The 31 audited declarations use only `propext`, `Classical.choice`, `Quot.sound`; no `sorryAx` or custom axiom. |
| Avoid mathematical smuggling | Analytic PDE existence, Sobolev estimates and continuum convergence are distinguished from abstract Lean algebra. Required source laws, feasible-class coverage, storage and blow-up returns remain explicit and OPEN where unpaid. Finite sampled diagnostics do not become uniform source or continuation theorems. |
| Fix experiment logic | Integrating-factor stages, aliasing, real spectra, restriction, supplied-source signs, capacity and localization scaling, residual return, causal frequency response, actual two-stage composition, stochastic conservation and terminal recording are repaired. Failed sweep cases cannot disappear into a positive overall verdict. |
| Verify repairs | The full suite passes all **110 tests**, including slow tests and independent analytic/false-target controls. Python compilation and `git diff --check` pass. |
| Check complete figure production | All **15** manifest producers have verified PNG/NPZ outputs at their default settings, combining the comprehensive run, corrected BU-29 rerun and final fast run. Both 18-case sweeps have complete coverage with zero failed configurations. All six fast NPZ outputs are byte-identical on rerun. See [artifact evidence](completion_artifacts.json) and the qualifications below. |
| Distinct adversarial pass | [Parent self-review](mathematics_self_review.md) checks signs, endpoint Sobolev exponents, semigroup/continuation applicability, zero cases, μ uniformity, discrete/continuum scope and diagnostic currencies. No independent agent review is claimed. |
| Keep the paper unchanged | No manuscript, figure manifest, historical figure, copied foundational paper, or separate unconditional campaign file is changed. |
| Stay within scope | No new source campaign, scout, consumer extension, paper expansion, Clay claim or later-paper import. Repairs remain uncommitted for review; nothing has been pushed. |

## Reproduction

The Python checks used the dependency environment
`/tmp/navier-scaffold-review-venv`, with `PYTHONPATH=python/src`,
`MPLCONFIGDIR=/tmp/navier-scaffold-mpl`, `OMP_NUM_THREADS=1` and
`OPENBLAS_NUM_THREADS=1`. Run from the repository root:

```sh
python3 -m pytest python/tests
python3 -m compileall -q python/src python/scripts python/tests
git diff --check
```

With the repository's Lean toolchain, run from `formal/`:

```sh
lake build
lake env lean ScaffoldAxiomAudit.lean
```

The environment versions were Python 3.12.3, NumPy 2.4.4, SciPy 1.17.1,
Matplotlib 3.11.2 and pytest 9.0.3. Use a Python environment containing
`python/requirements.txt`; the commands above refer to that environment.

The full default artifact command is:

```sh
python3 python/scripts/ns_paper_figures.py --include-slow --include-optional --out-root /tmp/navier-scaffold-completion-artifacts
```

During the review, the old BU-29 invocation was deliberately terminated
with SIGTERM after its duplicate rank implementation was found to need
repair. The comprehensive runner therefore reports that entry's exit -15
and exits 1; that receipt is retained and is not described as a green run.
All its other 14 entries passed. The full corrected BU-29 sweep was run
with:

```sh
python3 python/scripts/ns_paper_figures.py --include-slow --ids fig_nsbu29_antilocalization_sweep --out-root /tmp/navier-scaffold-corrected-sweep
```

The six final fast producers were also run with `--verify-determinism`
into `/tmp/navier-scaffold-final-fast`. Their two runs had identical NPZ
hashes. Final artifact evidence selects these fast outputs, the corrected
BU-29 output, and the other slow/optional outputs from the comprehensive
batch. The [JSON record](completion_artifacts.json) stores their exact
paths, hashes, configurations, coverage and current source hashes.
Historical manuscript artifacts have been preserved.

## Full artifact coverage

| Manifest producer | Class | Checked result |
| --- | --- | --- |
| NS-14 energy ledger | fast | PNG/NPZ present; deterministic NPZ; finite residual diagnostic |
| NS-15 UV completion | fast | PNG/NPZ present; deterministic NPZ; finite resolved/high-band comparison |
| NS-16 μ recovery | fast | PNG/NPZ present; deterministic NPZ; finite-grid low-mode error decreases in sampled μ values |
| NS-18 exponent harness | fast | PNG/NPZ present; deterministic NPZ; interpolation control distinguished; not the global-energy proof |
| BU-12 toy Zeno gallery | fast | PNG/NPZ present; deterministic NPZ; explicit geometric/harmonic examples |
| SPT legal/illegal control | fast | PNG/NPZ present; deterministic NPZ; generated baseline passes and localized control fails |
| BU-21 capacity numbers | slow | PNG/NPZ present; supplied-force convention and nested composition; unavailable initial ratios retained explicitly |
| BU-23 hypothesis screen | slow | PNG/NPZ present; sampled raw tolerance result and finite-prefix fit explicitly scoped |
| BU-28 localization | slow | PNG/NPZ present; measured localization/channel statistics only |
| BU-29 corrected sweep | slow | PNG/NPZ present; all 18 configurations; no failures; complete coverage; channel-normalized threshold passes |
| CHAN-00 rank | slow | PNG/NPZ present; all three seeds; no failures; pooled and windowed ranks separated |
| CHAN-01 delocalization | slow | PNG/NPZ present; all three seeds; no failures; full-snapshot bound includes discarded residual |
| SPT sweep | slow | PNG/NPZ present; all 18 configurations; no failures; complete coverage; raw-factor threshold passes |
| BU-04 response | optional | PNG/NPZ present; extracted-power sign named; negative as well as positive measured shells retained |
| BU-05 passive model fit | optional | PNG/NPZ present; actual optimizer succeeds; finite fitted model and error; no nonlinear bridge certificate |

## Artifact interpretation checked

The largest raw SPT concentration factor is approximately 5.78575, versus
1.15715 for BU-29's channel-normalized factor. Both have largest burst
ratio approximately 1.17130. These are finite sampled threshold results.
Common shells cover 18 cases; the highest listed shell covers the nine
N=48 cases and is empty after strict de-aliasing. Its zero rank and local
fraction express the absent-band convention, not a high-frequency theorem.

The channel-span bound at j=4 is approximately 0.002361. The worst actual
snapshot residual fraction is approximately 0.186124, and the returned
full-snapshot local-fraction bound is approximately 0.225868. Reporting the
span value as a full-field bound would discard a material term; the
repaired output keeps these distinct.

BU-21 has two intentional NaNs in `r_tilde`: the nonexistent initial
ratio and the ratio with a zero preceding capacity. They are unavailable
step factors, not numerical evolution failures. All other numeric fields
in the selected 15 artifacts are finite. BU-23 starts at the first positive
capacity; its sampled tolerance and 5% fit margin establish no infinite
raw summability or polynomial-growth law. The exact instantaneous nested
field composition is not a packaged-history source estimate.

The observed SGS response is mostly negative in two of the three tested
shells under the named extracted-power convention. BU-05 fits the remaining
eligible shell. Its weighted RMSE is approximately 0.001766 while the
unweighted RMSE is approximately 0.2444, and its real response visibly
misses much of the observed high-frequency data. Nonnegative coefficients
certify the chosen rational model's passivity, not the actual nonlinear NS
bridge. Choosing the smallest computed weighted error among four fits is
not a global optimum or a uniform constitutive law.

Representative generated plots were inspected as well as their numeric
arrays. These interpretation limits are part of the review, even though
all selected producers successfully wrote their required artifacts.

## Mathematical limits of completion

Completion of this review does not close the original-NS supplier or
blow-up-to-frontier obligations. It repairs and verifies the existing
conditional scaffold. The [hyperviscous proof](proof_skeleton_ns_alpha.md)
is for fixed μ>0; [convergence](smooth_interval_recovery.md) assumes a
smooth original NS reference through the supplied endpoint. The
[scalar-symbol result](passive_low_frequency_closure.md) does not derive
its constitutive premises from memory autonomy or fixed-mode passivity.
The unchanged paper requires the specific scope/convention/figure updates
listed in the review during its later editing phase.
