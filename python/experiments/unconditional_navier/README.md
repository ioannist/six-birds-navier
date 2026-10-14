# Reproducible construction controls

Run `python3 python/experiments/unconditional_navier/spectral_controls.py
--output python/artifacts/unconditional-navier/controls.json` from the root.
Python's standard library suffices; Fraction arithmetic gives exact finite
identities without floating-point evidence claims.

The domain is a finite periodic transverse Fourier control carrier. These
are deliberately separate from the original R3 Navier trajectory. The
parent's analytic annulus and half-open-path constructions supply precise
logical controls; the finite computation only checks sample identities.
New equation experiments must state their original-domain/limit bridge,
signed interaction term, cutoff errors and source applicability before
being promoted as a leverage rung.

Run `python3 python/experiments/unconditional_navier/source_scaling_controls.py`
for exact critical-source profile diagnostics. For
`u=tau^(-beta) U(x/tau^alpha)`, rational exponents compare energy,
enstrophy and four conditional source doses. The parabolic profile has
logarithmically divergent critical doses despite bounded energy and
integrable enstrophy. The historical energy profile has divergent critical
doses as well. These are profile scaling controls, not asserted solutions
of the unforced PDE. No Navier criterion is falsified by a non-solution.

Run `python3 python/experiments/unconditional_navier/heat_ledger_controls.py`
for exact ORIGINAL Schwartz initial-state heat-transfer controls. Rational
polynomial/Gaussian arithmetic checks solenoidality, the heat generator,
the all-scale signed transfer formula and its positive factor. A
fixed-energy spatial rescaling gives transfer growth R^(5/2), versus R^2
for the viscous complement. This rejects the specified instantaneous
energy/viscosity bound even at fixed viscosity; it supplies neither a
single-data singular trajectory nor a terminal-uniform signed source.
The analytic bridge and exact scope are recorded in the project's
research notes, which are kept outside this repository.

Run `python3 python/experiments/unconditional_navier/heat_window_controls.py`
for exact directional-derivative and source-price controls. Independent
rational polynomial and direct Gaussian calculations give a negative
viscous cubic derivative at every positive heat scale. The written
small-amplitude original-data bridge shows why a physical derivative cannot
be omitted; that experiment leaves the quartic pressure coefficient uncomputed.
The separate round09 one-receiver experiment below computes it at beta=3/8.
Exact lag-band charge exponents show that the native absolute energy price
is summable in the tested algebraic parabolic pattern. That pattern is
not asserted to be a realized Navier trajectory or disjoint source sequence.

Run `python3 python/experiments/unconditional_navier/source_window_controls.py`
for exact native interpolation exponents, reflected triad geometry and the
correct normalized split-minimizer diagnostic. These rational checks do
not construct source events, a cadence or a signed PDE source bound.

Run `python3 python/experiments/unconditional_navier/vorticity_controls.py`
for exact original Gaussian transport/stretching branches, their signed
Fourier conversion, and actual swirl/strain-jet applicability controls.
They supply no time-integrated source or general symmetry applicability.

Run `python3 python/experiments/unconditional_navier/geometry_controls.py`
for exact original directional and bending jets. Plain Fraction polynomial
arithmetic independently compares direct curl-PDE material prolongation
with the complete retained bending balance. A separate compactifiable
solenoidal shear jet has POSITIVE complete pointwise viscous bending source;
this refutes pointwise negativity, not global integral negativity or
regularity. The small energy/enstrophy-product threshold is checked as an
exact rational constant; its mathematical class theorem is in round08.

Run `python3 python/experiments/unconditional_navier/quartic_controls.py`
for the exact one-heated-receiver quartic at beta=3/8 with full Leray
projection, physical/Fourier raw cross-checks, beta0 cancellation and a
rational radical sign certificate. It gives a negative INITIAL coefficient;
it does not sign the full trajectory source. Run
`python3 python/experiments/unconditional_navier/bending_integrated_controls.py`
for four actual solenoidal Schwartz Gaussian integrated-viscosity controls
and the p=1,2,4 compactifiable shear normalization controls. These two
helpers require SymPy; the recorded runs used version1.14.0. Exact artifacts
are written under python/artifacts/unconditional-navier/. The four negative
Gaussian values do not establish a universal integrated sign. The separate
round10 analytic finite-scale curl/heat construction proves a positive
integrated coefficient on another actual Schwartz state; no floating
computation or auxiliary heat path is promoted to a Navier trajectory.

Run `python3 python/experiments/unconditional_navier/critical_pressure_source_controls.py`
for the exact finite Fourier pressure triad, positive boosted work coefficient,
once-only envelope and fixed-input mass scaling checks. This requires SymPy
and writes `python/artifacts/unconditional-navier/critical-pressure-source-controls.json`.
It checks finite algebra only. Whole-space curl compactification, the slow
pressure mode, both nonlocal pressure crosses and actual local NS application
are proved separately in the round67 originals and parent audit. The resulting
fixed-E/fixed-I control refutes only the specified instantaneous norm-only
source bound; it does not construct a blowup or a positive regularity estimate.
