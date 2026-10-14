# Conditional recovery on a smooth Navier Stokes interval

Fix a periodic box, ν>0, α≥5/4, common smooth solenoidal initial data,
and common smooth forcing. Let v be the classical Navier–Stokes solution
on a **supplied smooth interval [0,T]**, and let u_μ be the hyperviscous
solution with coefficient μ>0. The fixed-μ existence theorem is proved
in [the hyperviscous dossier](proof_skeleton_ns_alpha.md).

For every integer m≥3 there are constants K_m and μ_m>0, depending on
the reference v, T, the box and m,α,ν, but not on μ, such that

    sup_{0≤t≤T} ||u_μ(t)-v(t)||_{H^m} ≤ K_m μ,    0<μ≤μ_m.

Thus u_μ→v in C([0,T];H^{s'}) for every finite s'. This restores the
existing conditional recovery claim. Smoothness of v through T is a
premise, not a consequence of hyperviscosity or the convergence estimate.
The claim provides no original-NS continuation beyond T.

## Actual difference equation

Put w=u_μ-v, with w(0)=0. Subtract the two projected equations and write
the regularization as dissipation on w plus a source evaluated on v:

    ∂t w + νΛ²w + μΛ^{2α}w
      + B(v,w)+B(w,v)+B(w,w) = -μΛ^{2α}v.

The common forcing cancels. If the data/force have nonzero spatial mean,
both solutions have the same mean, and the translating-frame adapter in
the existence dossier applies to both. Its spatial Sobolev norms are
translation invariant. No high Sobolev bound for u_μ, uniform in μ, is
assumed.

## Difference energy estimate

For an integer m≥3, ordinary derivative Leibniz estimates, transport
cancellation and H^m⊂W^{1,∞} on the fixed three-dimensional torus give

    |<B(v,w)+B(w,v)+B(w,w),w>_{H^m}|
      ≤ C_m (||v||_{H^{m+1}}+||w||_{H^m}) ||w||_{H^m}².

At integer order the energy uses the sum of squared L² norms of derivatives
through order m, an equivalent standard H^m norm. Its fixed equivalence
constants are included in C_m and the returned K_m.

Here is the derivative bookkeeping. In B(v,w), the highest derivative
transport v·∇D^βw cancels against D^βw by div(v)=0. The remaining terms
are bounded by the commutator/product estimate
||∇v||_∞||w||_{H^m}+||v||_{H^m}||∇w||_∞, paired with ||w||_{H^m}.
In B(w,v), the product bound uses ||w||_{H^m}||v||_{H^{m+1}}.
For B(w,w), its highest transport cancels and the remainder is bounded
by C_m||∇w||_∞||w||_{H^m}². The embedding for m≥3 pays all displayed
supremum norms. The Leray projector drops from the pairing with the
solenoidal w. The estimate retains all three nonlinear difference terms.

Define X=||w||_{H^m}², M=sup_{[0,T]}||v||_{H^{m+1}} and
F(t)=||Λ^{2α}v(t)||_{H^m}. Pairing the equation gives

    (1/2)X' + ν||Λw||_{H^m}² + μ||Λ^αw||_{H^m}²
      ≤ C_m(M+sqrt(X))X + μF(t)sqrt(X).

The last term is bounded by (X+μ²F(t)²)/2. In the region X≤1, therefore,

    X' ≤ A_m X + μ²F(t)²,   A_m=1+2C_m(M+1).

The constants in the product estimate are geometric/Sobolev constants;
the coefficient A_m contains only the fixed smooth reference and the
bootstrap threshold, not an independently supplied bound on u_μ.

## Bootstrap and rate

Set

    K_m² = exp(A_m T) ∫_0^T F(t)² dt,
    μ_m = 1/(2 max(1,K_m)).

Up to any first exit from X≤1, Grönwall and X(0)=0 imply

    X(t) ≤ μ² ∫_0^t exp(A_m(t-s))F(s)² ds ≤ μ²K_m² ≤ 1/4.

Continuity rules out a first exit at X=1. Since u_μ exists on the full
interval for each μ>0 and the reference is smooth there, the estimate
holds through T. Taking square roots proves the rate, including K_m=0.
Sobolev inclusion/interpolation returns every finite fractional s' from
an integer m≥max(3,s'). All constants are uniform as μ decreases through
the stated small interval.

The displayed O(μ) rate uses reference regularity H^{m+2α}, as well as
H^{m+1}. Smooth data and a smooth reference pay this at every finite m.
If only finite reference regularity is supplied, this derivative cost
must be retained; an O(μ) claim at arbitrary higher order cannot be made
by hiding those derivatives inside “standard stability.”

## Mechanization and experiment boundaries

`Fluids.Difference` proves abstract Hilbert norm differentiation and
dissipative forcing algebra. It does not encode the Sobolev commutator
estimate or the continuum PDE existence/continuation argument above.
The continuum recovery proof is analytic, not a Lean PDE theorem.
Fixed-grid μ sweeps test a different, finite-dimensional consistency
question and do not supply its smooth-reference premise.

An exact shear mode u=A exp(-(ν|k|²+μ|k|^{2α})t)e sin(k·x), e·k=0,
has zero advection. Its completed/reference difference divided by μ
tends to A|k|^{2α}t exp(-ν|k|²t)e sin(k·x). This checks the sign and
first-order μ rate independently of a nonlinear numerical fit.
