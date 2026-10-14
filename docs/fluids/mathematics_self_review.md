# Adversarial self review of the Navier Stokes corrections

The repaired abstract scaffold and supporting derivations retain their
existing scope. This separate pass checks the mathematical changes against
the original statements and return obligations. It is parent self-review,
not independent external or agent review.

## Domains and quantifiers

The paper is a periodic discrete/Hilbert-space decision scaffold. A Fourier
ball on R³ is never used as a finite-dimensional object. The sharp shell
choice is explicit; smooth overlapping filters cannot silently inherit the
projection return. No stronger later Navier result or unconditional closure
premise is inserted into the old scaffold.

The tame-factor adapter uses C(n)≥0 at every step, already part of the
paper's setup. The old stronger direct-recurrence exports are preserved.
The new Parseval proof expands the actual Gram matrix; it does not accept
field-energy equality as an unexplained hypothesis. Basis feasibility is
supplied explicitly for the restricted ICAP theorem. Empty spatial/channel
sets give zero finite sums; no nonempty assumption or division is hidden.
At fixed channel delocalization constant, m_j 2^{-3j}→0 is needed to
return η_j→0. Merely having finitely many channels at each shell is
insufficient. Exclusion of tests with a fixed positive local energy
fraction is quantified at sufficiently large shells; it does not pay a
continuum gain representation or pressure/strain limit.

The original No-Zeno theorem only needs unbounded partial sums and an
increment inequality. Monotonicity is a separate result requiring
nonnegative increments. Its dyadic wrapper indexes crossing times by
j and samples the increment certificate at 2^j; it does not infer dyadic
crossing hypotheses from a full sequence without an adapter.

## Source signs and exact counterexample

With NL=P div(u⊗u), the actual equation has RHS -NL. The supplied closure
is NL(v)-NL(u) after low projection. BU-21 therefore takes a difference
of the already negative numerical nonlinear function, whereas BU-04/05
use its opposite sign to fit extracted power. Positive supplied work and
positive-real extracted response must not be conflated.

For the triad k=(4,0,0), l=(-3,1,0), q=(1,1,0), the stated polarizations
are solenoidal. The projected sum interaction is +(3/2)c sin(q·x) in NL,
so its supplied closure is negative. Pairing with -εc sin(q·x) makes
positive power +(3/2)ε Vol. Since |c|²=2 and the spatial mean of sin² is
1/2, shell norm squared is ε² Vol. The ratio 3/(2ε) is correct. Mixed
low/high modes miss the tested shell. The datum is a genuine smooth
periodic velocity, not a formal unattainable source assignment.

The repaired shell estimate counts finite output frequencies before applying
Leray. It does not assert an L¹ Leray bound. It controls absolute power by
native energy but does not invent a shell-energy lower bound. No refutation
of all nonlinear functionals or temporal capacity bounds follows.

## Hyperviscous endpoint and forcing

At α=5/4, H^{1/4} embeds into L^{12/5}, and H^{5/4} into L^{12} in three
dimensions. The two derivative-product interpolation weights sum to one;
the product lies in L^{12/11}, dual to L^{12}. Thus higher nonlinear work
is bounded by a(t) X_m^{1/2}D_m. Young leaves a(t)²X_m/μ, whose integral
is paid by native dissipation. No critical L∞ embedding is assumed.
The local Duhamel map loses one derivative in advection and regains it
with kernel (μt)^(-1/(2α)); its exponent is below one. Its restart time
depends on the bounded H^m ball at fixed μ, exactly what the continuation
argument supplies. No lifespan uniform as μ→0 is inferred.

For α>5/4, the fixed torus's positive smallest frequency controls the
lower dissipative order, with a domain/α constant. Forcing is controlled
in Λ^{-α} at the base energy level and in H^m at the higher level. Smooth
forcing with nonzero mean returns through the explicitly known translating
frame; its gradient part is put in pressure. The supremum plus full
base-dissipation bound has the required factor two. Constants deteriorate
with μ→0; the proof does not assert an original NS theorem at μ=0.

## Crossing and divergence returns

A finite continuous storage crossing attains the threshold on a compact
interval. Passivity and monotonicity of positive work then give WORK.
Infimum notation alone, or an absolute shell threshold already exceeded
at the initial crossing, gives no fresh-work quantum. An infinite frontier
is allowed in the definition until No-Zeno is proved.

A capacity upper bound gives a lower reciprocal bound. A convergent such
lower bound decides nothing about the actual reciprocal series. The
constant-capacity control checks this quantifier. The conditional unforced
Clay-facing chain now uses native energy to supply a uniform shell budget,
with zero data handled separately. All PDE-hard source and storage-return
hypotheses stay explicit. Fixed-time H^k tails vanish with cutoff; the
old opposite tail limit is removed.

## Numerical evidence limits

Integrating-factor stages reduce to RK4 when L=0 and to the exact linear
semigroup when nonlinear/forcing terms vanish. Forced-linear and quadratic
analytic controls discriminate the repaired fourth-order scheme from the
old wrong weights. Strict de-aliasing removes the 3K=N collision. Nested
Fourier restriction drops coarse Nyquist planes without merging one side
of a conjugate pair; same-grid restriction remains the identity.

Relative rank and capacity do not depend on an absolute amplitude floor.
The certificate fails on missing required shell ranks. Pooled rank remains
an energy statistic; no exact feasible subspace is inferred. The local
triangle inequality charges a concentrated orthogonal residual and returns
the measured span bound only with that charge.

The Debye sign is checked against a stable ODE, rather than generating
both test data and predictions with the same implementation. Signed-frequency
data use absolute frequencies for positive-rate initialization; an empty
memory sum returns the constant response at every requested frequency.
Empirical
closure-fit improvement is not a conserved quantity. Empty high-frequency
bands give an unavailable result rather than a forced positive slope.
Floating-point tests and fits remain diagnostics, not interval certificates
or continuum proofs.

## Conditional recovery and scalar-symbol recheck

The difference equation keeps both mixed transport terms and the quadratic
term in the difference. At integer m≥3, H^m embeds into W^{1,∞} in dimension
three; the reference's H^{m+1} norm pays the derivative in B(w,v). Testing
the regularization on the reference costs H^{m+2α}. The source sign is
negative, but its absolute energy contribution is μF sqrt(X). Young and
the first-exit argument close X≤μ²K²≤1/4 without a uniform higher bound
on the completed solution. The zero-source and T=0 cases also return.
Every returned fractional Sobolev order is supplied by a larger integer
order. This is a conditional analytic proof, not an original-NS existence
theorem or a Lean PDE limit theorem.

The low-frequency symbol conclusion uses positivity for all sufficiently
small positive r. The polynomial r²-r/2 checks that positivity only at
realized torus modes is insufficient, even for a smooth scalar symbol.
The Lean counterexample proves its derivative and mode-region positivity.
Taylor's remainder requires the separately stated C² bound. A zero
coefficient is permitted; the total existing-plus-closure coefficient is
ν+η. Memory autonomy and scalar isotropy are not derived from role labels.

## Completion-audit experiment controls

Terminal-time tests use exact shear solutions, including horizons shorter
than one nominal step. Both dimensions return the actual requested terminal
time and its state. Hermitian input checks prevent silently discarding an
imaginary component, and nonfinite states cannot become successful runs.
The de-aliasing-off control removes only ambiguous Nyquist faces; it is
still a control with aliasing, not an energy-conserving spectral method.
The dyadic toy's one-mode equation has an independent exponential solution;
it checks short/fractional horizons and sparse output recording. Its signed
damping coefficients remain allowed for the intended active controls. All
five default toy configurations reach their requested t=0.3 endpoint.

The two-stage closure returns the actual high-to-middle response before
restriction to the coarse grid. Its independent nonlinear example makes
the formerly omitted response nonzero. Ordered checkerboard particle
updates use disjoint bonds on even periodic lattices; exhaustive four-site
states and deterministic probabilities verify species conservation. The
signed-current interpretation is not inferred from an unsigned phase rate.

Constant raw mismatch with decaying depth-normalized mismatch is a
false-target control for the raw summability screen. Empty capacities do
not yield a measured step factor. Every finite prefix admits both summable
and nonsummable continuations, and a fixed relative per-step margin incurs
an exponential accumulated loss. Diagnostic PASS/PLAUSIBLE labels retain
their finite sampled meaning.

Failed-configuration controls exercise both sweep aggregators: favorable
surviving cases cannot produce an overall positive verdict when coverage
is incomplete. The SPT sweep screens concentration and burst; its m95 is
reported as a diagnostic, whereas the separate certificate tool explicitly
checks m95 when `require_m95` is supplied. BU-29's channel-normalized factor
is labeled separately from the raw concentration factor used by that
certificate. No finite sweep supplies a uniform actual-history theorem.
At N=48, strict de-aliasing retains component frequencies at most 15,
so the reported shell j=5 (|k|≥32) is empty. Its zero concentration/rank
is an absent-band convention, not evidence of high-frequency
anti-localization. The corresponding finite-grid qualification applies
to both default sweeps.

The full current Python suite passes all 110 tests. All 31 original-plus-
adapter Lean declarations have only the standard `propext`,
`Classical.choice` and `Quot.sound` dependencies. Semantic returns remain
distinct from the successful build and transitive-axiom checks.

## Review outcome

The corrections do not materially downgrade any main abstract theorem.
The physical results repaired here are the classical fixed-μ>0
hyperviscous theorem at α≥5/4 and its convergence to a supplied smooth
original reference interval; both are analytic, not mechanized. The rejected
SGS estimate and unsupported continuum no-go shortcut are confined to
their actual scope. The manuscript needs the recorded wording and artifact
updates in its separate editing phase. No new construction campaign,
consumer extension, paper expansion or unconditional target gain is claimed.
