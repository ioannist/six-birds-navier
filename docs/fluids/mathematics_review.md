# Mathematics review of the Navier Stokes scaffold

The abstract results supporting *Six Birds for Navier–Stokes: A Mechanized
No-Zeno Scaffold* are correct at their stated discrete or Hilbert-space
scope. No material downgrade to the title or main abstract theorems is
needed. The supporting mathematics and numerical implementations required
the corrections below. The original manuscript remains unchanged for its
later editing phase.

This review covers the original `paper/main.tex` section chain, all nine
original `formal/Fluids/*.lean` modules, their exported declarations and
transitive axioms, the associated `docs/fluids` derivations, and the numerical
operators, solvers, capacity, localization, rank and frequency-response
code used by the scaffold. The later sibling Navier theorem, historical
implementation copies, and the separate unconditional research campaign
are outside this review. No stronger conclusion is imported from them.
The pre-review state is preserved in commit `b2c2316`; rewritten derivations
also have local original copies under the ignored `review/` directory.

## Exact theorem coverage

| Component | Strongest verified conclusion and its premises |
| --- | --- |
| Energy | Assuming quadratic inner-product cancellation, PSD forms and nonnegative coefficients, the abstract RHS has nonpositive pairing with the state. This is not a formal PDE energy identity. |
| Multiplier | A function vanishing at zero, differentiable there, and nonnegative in a right neighborhood has nonnegative derivative. No NS multiplier expansion is derived from this algebra. |
| Difference | Norm-square differentiation and dissipative forcing inequalities hold along a supplied differentiable curve satisfying the abstract RHS relation. The analytic continuum recovery theorem now has a separate explicit Sobolev bootstrap; it is not mechanized by this algebra. |
| No-Zeno | Increment lower bounds and unbounded partial sums imply that crossing times have no finite upper bound. Nonnegative increments additionally imply monotonicity. |
| Work and capacity | Positive capacity and the supplied work/cap inequality permit division; unbounded ratio sums then give No-Zeno. The storage construction and continuous time integrals are separate analytic inputs. |
| Route capacity | The polynomial recurrence gives the stated finite-sum bound. Bounded nonnegative mismatch sums give its uniform polynomial version. The actual tame-factor adapter needs the nonnegative capacity sequence stipulated in the paper's setup. |
| Aggregation | Positive parts are subadditive under finite sums, and per-atom bounds aggregate with the number of atoms. A continuum NS atom decomposition is not constructed. |
| Anti-localization | Per-channel cell bounds give the coefficient-energy inequality. The orthonormal Gram identity gives Parseval and the return to actual field energy. Uniformity over all cells requires the channel bound for each cell. |
| ICAP no-go | A discrete diagonal quadratic gain certificate on a feasible class containing every basis vector bounds each positive diagonal gain. This is not an automatic theorem for the nonlinear packaged NS history map. |

The missing Lean connections have been supplied without changing the old
exports: `route_capacity_poly_bound_of_tame_factor` and its bounded-sum
version; `synthesis_energy_eq_coeffEnergy_of_orthonormal` and its cell-mass
corollary; and `icap_implies_pointwise_pos_bound_on_feasible_basis_tests`.
The derivative-operator return and fixed-box passivity counterexample are
also proved in `Difference` and `Multiplier`, respectively.
They prove the adapters from structural hypotheses, rather than asking
for their desired conclusions as proof-carrying inputs.

## Mathematical corrections

1. **The bridge must use actual projections.** Smooth overlapping
   Littlewood–Paley multipliers are not idempotent. Defining their ranges
   by fixed points and then setting the high state to `(I-G_j)u` did not
   give the asserted actual-path return. `ns_bridge_object.md` now fixes
   the sharp orthogonal shells already used in the scaffold diagnostics,
   proves the high-range invariance and local strong uniqueness needed
   for causality, and returns the readout on an actual strong path.
   The input is the whole resolved history, not just the shell port.

2. **Supplied and extracted power had opposite conventions.** The actual
   resolved forcing is `G_j NL(G_j u)-G_j NL(u)`. Since the numerical
   nonlinear function already returns `-NL`, BU-21 must take its direct
   incremental difference. The extra minus has been removed. BU-04/05
   retain an explicitly named extracted-power/LHS diagnostic convention;
   fitted passivity there is not positive supplied-work ICAP.

3. **The full SGS shell-only strain estimate is false.** The high–high
   source can populate a shell even when that shell is initially absent.
   `cap_bound_attempt.md` gives a smooth solenoidal three-mode datum with
   bounded full gradient and positive supplied power proportional to ε,
   while shell squared norm is proportional to ε². Its exact ratio is
   3/(2ε). The claimed uniform bound by `||∇u||∞ ||u_j||₂²` is therefore
   withdrawn. The correct finite-shell estimate is
   `||y_j||₂ ≤ C 2^{5j/2} ||u||₂²`, giving an absolute energy-based power
   bound, with no hidden lower bound on shell energy. This is a scoped
   counterexample to that estimate, not a no-go for all functionals or
   cumulative history bounds.

4. **Capacity comparisons had reversed quantifiers.** A polynomial upper
   bound with exponent at most one guarantees reciprocal divergence.
   An upper bound with a larger exponent or exponential growth does not
   establish convergence of the actual reciprocal series. Constant
   capacity is an immediate counterexample. Failure of the comparison
   certificate proves neither Zeno nor impossibility of a sharper CAP.
   The affected capacity, divergence and work-quantum notes are corrected.

5. **Crossing infima and PDE tails require their own returns.** Continuous
   storage and finite nonempty attained crossings justify the work quantum;
   merely defining an infimum does not. The frontier may be infinite and
   must not be defined using a maximum assumed to exist. At every fixed
   smooth time an H^k tail tends to zero with the cutoff, rather than to
   infinity. The actual blow-up-to-storage-frontier implication remains
   open. Native unforced energy supplies the uniform shell budget needed
   to turn the conditional linear ICAP bound into linear total capacity.
   Proposed blow-up-to-passivity-failure and growing-mismatch implications
   are explicitly OPEN as well; diagnostics do not prove those returns.

6. **The hyperviscous endpoint proof is repaired at its original strength.**
   The old high-norm-only interpolation produced a superlinear scalar
   inequality and did not close by Grönwall. The corrected proof in
   `proof_skeleton_ns_alpha.md` uses the integrable native-energy coefficient
   `||Λ^{5/4}u||₂²`, finite-exponent periodic Sobolev embeddings and derivative
   product interpolation. It retains α≥5/4 for fixed μ>0 and ν≥0, covers
   smooth forcing with an explicit mean/gradient adapter, and supplies the
   higher-energy estimate needed for continuation. It does not provide
   bounds uniform as μ→0 or a theorem for μ=0. The supporting primary
   source is [Tao's critical energy method, Section 2](https://msp.org/apde/2009/2-3/apde-v2-n3-p04-p.pdf);
   the fixed-torus and forcing derivation is written out locally.

## Additional conditional results restored

The original fluids anchor's Theorem 3 now has a full proof in
[smooth_interval_recovery.md](smooth_interval_recovery.md). With common data
and forcing and a supplied smooth original NS reference on [0,T], the
regularized solutions converge in every finite Sobolev order with an explicit
O(μ) rate for sufficiently small μ. The completed solution's high norm is
controlled by a closed bootstrap, not assumed uniformly bounded. The rate's
reference derivative cost H^{m+2α} is visible. This preserves the intended
conditional result and does not prove original NS continuation.

[passive_low_frequency_closure.md](passive_low_frequency_closure.md) supplies
the scalar time-local symbol/Taylor implication at its actual applicability
scope. Passivity at a fixed torus's Fourier modes cannot replace positivity
near zero; a smooth positive-mode symbol with derivative -1/2 at zero is now
also a Lean counterexample. A nonnegative viscosity coefficient can vanish,
and autonomous memory lifting does not imply a velocity-only time-local law.
The correct total leading coefficient in an equation already containing νA
is ν+ℓ'(0). The historical anchor plan is retained with an explicit correction
notice pointing to these current derivations.

## Numerical corrections and controls

Both spectral solvers had incorrect nonlinear weights in their time stages.
They now apply classical RK4 to the integrating-factor variable; the 2D
historical `step_etdrk4` name is retained with an accurate docstring. Exact
forced-linear and quadratic scalar solutions independently check fourth-order
convergence. The 2D solver also masks its supplied initial data and forcing
before using them in stages.

The de-aliasing cutoff is strict: a boundary with 3K=N lets two retained
modes alias to a retained mode. Tests include grids divisible by three and
check the semi-discrete nonlinear energy cancellation. Fourier restriction
now drops ambiguous coarse Nyquist faces, preserving Hermitian symmetry,
real fields, Fourier amplitude normalization and nested restriction.

The sampled work/energy ratio no longer adds an absolute denominator epsilon
that understated capacities and broke amplitude invariance. It is labeled
as a ratio on one sampled window, not an extremum over feasible histories.
Local fractions and effective ranks are rescaled to avoid losing tiny
nonzero fields. Required rank checks reject missing, negative or fractional
shell ranks, and invalid concentration/baseline ranges are rejected.

Pooled m95 measures uncentered sampled squared energy, not exact feasibility
or centered statistical variance. A small discarded component can be highly
localized. Channel diagnostics now label eta_hat as a bound on the retained
span, record actual per-snapshot residual energy, and return a full-snapshot
bound using the local L² triangle inequality. A concentrated orthogonal
residual is the false-target control; fixed 95% retention cannot prove
asymptotic anti-localization.
A finite channel count in each shell also does not imply asymptotic
anti-localization without the required sub-volume growth condition.

The causal Debye response now uses `1/(1+iω/b)` for NumPy's negative-exponent
Fourier transform. An independent stable first-order ODE transfer law checks
the phase and the fit, including negative frequencies and the constant
zero-memory case. Positive-real fit coefficients certify the chosen
rational model only, not the actual nonlinear NS bridge.

The slow closure test incorrectly required averaging to improve a particular
empirical fit. After the mathematically necessary de-aliasing fix, the old
high band has no retained modes on N=12. The test now verifies finite
numerical output and the correct unavailable-band result; an exact Laplacian
symbol separately tests the fit and slope calculations. A slow test's
repository path was also corrected. No solver parameter was tuned to recover
an unsupported favorable empirical conclusion.

## Additional experiment logic repaired

- The percentile frontier now counts only energy strictly above the retained
  shell. Parseval normalization gets its spatial axes from the shell mask,
  avoiding a scalar N=2 or N=3 grid being mistaken for a vector field.
- The ordered checkerboard Markov update requires an even periodic lattice;
  on an odd lattice overlapping bonds could destroy a particle. Exhaustive
  four-site controls check both species for every state and deterministic
  probability choice. Probabilities and initialization amplitudes retain
  their stochastic domain. The unsigned phase-rate proxy is not reported
  as a signed drift or as proof of a composed-kernel cycle affinity.
- Simulations end at the requested time, with a final short step when needed,
  and record the terminal state. Invalid nonreal Fourier spectra and nonfinite
  states/forces cannot silently become a passing real-field certificate or
  successful run. Aliased controls still remove ambiguous Nyquist faces.
  The supporting dyadic toy solver also returns the requested terminal time
  and state, and rejects a nonfinite energy/power diagnostic.
- Raw and depth-normalized mismatch trends have separate status. A decaying
  normalized series does not pay raw summability; finite sampled trends and
  a fixed 5% step tolerance do not prove an infinite polynomial bound. All-zero
  capacity data produce no fitted step factor and are INCONCLUSIVE.
- The original two-stage closure diagnostic omitted the high-to-middle
  response. It now retains both stages and checks exact nested composition.
- An ineligible passive fit writes the actual negative observations and an
  unavailable-fit status instead of omitting its artifacts or manufacturing
  a favorable model. Failed sweep configurations leave coverage incomplete
  and the overall result inconclusive. The sweep's legacy channel-normalized
  concentration factor is distinguished from the raw SPT certificate factor.
- Tiny nonzero fields can be energy-matched without an absolute amplitude floor;
  a zero field cannot meet a positive target energy. A PSD constructor rejects
  a negative scale rather than silently constructing an anti-passive operator.

## Changes reserved for the paper phase

The original manuscript files and figures have not been edited. Its later
revision should reflect the corrected supporting mathematics:

- Replace “summability implies no finite upper bound” in the introduction
  with **unbounded partial sums** or **divergence**.
- Identify the bridge shells as sharp projections and use the actual supplied
  force sign consistently. A smooth-filter variant needs a separate adapter.
- Limit the proved ICAP no-go to its finite diagonal model; any continuum
  wavepacket-to-strain passage needs an actual feasible test, gain and limit
  argument. Excluding those tests blocks that implication, but supplies no
  ICAP, WORK, DIV or continuation theorem by itself.
- Give the zero-shell convention for localization fractions and distinguish
  grid cubes, continuum balls, uncentered energy ranks and discarded residuals.
- Identify each formalization boundary accurately: neither storage hitting
  nor PDE continuation nor a continuum μ→0 limit is mechanized here; the restored limit proof is analytic.
- Regenerate numerical artifacts from the repaired solver and sign conventions
  when the paper is revised. Historical images have been preserved.

These are corrections to scope, conventions and diagnostics. The paper's
main abstract theorem statements and explicit non-claim of Clay regularity
remain intact. The unconditional campaign remains OPEN and unchanged.

## Validation

- `lake build` succeeds on Lean 4.8.0 with mathlib v4.8.0.
- `lake env lean ScaffoldAxiomAudit.lean` audits all 31 scaffold declarations:
  only `propext`, `Classical.choice` and `Quot.sound` occur transitively.
  There are no `sorryAx` or added custom axioms.
- The full Python suite, including slow tests, passes: **110 passed in 47.94s**.
- A small channel-diagnostic run exercises the recorded residual and full-field
  return after its repair. Python source compilation and `git diff --check` pass.
- A distinct [adversarial self-review](mathematics_self_review.md) checks the
  repairs' signs, endpoints, losses, quantifiers and semantic scope. This is
  self-review, not an independent agent review.

This report has been updated after the completion audit found previously
uncovered experiment logic. Its original first report and self-review remain
immutable in the private preservation tree. The completed
[requirement audit](completion_audit.md) records all 15 validated artifact
producers, both complete 18-case sweeps, six fast determinism checks,
source/artifact hashes and the exact diagnostic qualifications. No running
producer or failed stale invocation is counted as a successful validation.
