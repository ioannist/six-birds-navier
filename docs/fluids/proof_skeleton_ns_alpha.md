# Hyperviscous Navier Stokes energy proof

For fixed ν≥0, **μ>0** and α≥5/4, the periodic equation

    ∂t u + P((u·∇)u) + νΛ²u + μΛ^{2α}u = f,   Λ=(-Δ)^{1/2},

has a unique global smooth velocity for smooth divergence-free initial data
and smooth forcing. The argument below proves the a priori estimates used
with local strong existence and continuation. It is an analytic proof,
not a PDE formalization in Lean. Constants need not be uniform as μ→0.
In particular it supplies no regularity theorem for μ=0 in three dimensions.

## Periodic setting and forcing

Start on the fixed torus of side length L with zero-mean divergence-free u
and zero-mean solenoidal f. The Leray projector P and Λ commute with
spatial derivatives; P is self-adjoint. All Sobolev constants may depend
on L, the derivative order and α. Poincaré's inequality controls lower
orders by homogeneous norms on this zero-mean space.

General smooth forcing does not preserve zero mean. To cover it, first
project the force (its gradient part is absorbed into pressure). Let
m(t)=mean(u0)+∫_0^t mean(f(s)) ds and X'(t)=m(t), X(0)=0. Then
v(t,x)=u(t,x+X(t))-m(t) satisfies the same projected equation with forcing
P f(t,x+X(t))-m'(t), which is solenoidal and zero mean. Smoothness and all
finite-horizon spatial estimates return by reversing this explicit change
of variables. The pressure is reconstructed by the periodic Poisson
equation and a fixed zero-mean gauge. No mean-zero invariant is assumed
for an arbitrary force.

## Local strong existence and the continuation input

For integer m≥3, the periodic product estimate gives a locally Lipschitz
map B(u,u)=P div(u⊗u) from H^m to H^{m-1}. The linear semigroup
S(t)=exp(-t(νΛ²+μΛ^{2α})) obeys

    ||S(t)g||_{H^m} ≤ C(1+(μt)^{-1/(2α)}) ||g||_{H^{m-1}}.

Indeed, the nonzero-mode multiplier is bounded by
sup_{r>0} r exp(-μt r^{2α}) ≤ C_α(μt)^{-1/(2α)};
the inhomogeneous norm adds the displayed constant. Since 1/(2α)<1,
this kernel is integrable at t=0. The Duhamel map

    u(t)=S(t)u0+∫_0^t S(t-s)(f(s)-B(u(s),u(s))) ds

is a contraction on a ball in C([0,δ];H^m) for sufficiently small δ.
For fixed μ>0, its lifespan is bounded below in terms of the H^m ball,
the finite-horizon forcing norm and the fixed parameters. The same map
gives strong uniqueness and restarts at any preterminal time. Thus a
uniform preterminal H^m bound excludes a finite maximal lifespan by
restarting close enough to its endpoint. Applying the construction at
higher orders for smooth data, with uniqueness on overlaps, justifies
the smooth energy calculations below. No μ-uniform local lifespan is
needed or asserted.

## Native energy pays the critical coefficient

Nonlinear cancellation gives

    (1/2) d/dt ||u||₂² + ν||Λu||₂² + μ||Λ^αu||₂² = <f,u>.

Young's inequality bounds the forcing term by

    <f,u> ≤ (μ/2)||Λ^αu||₂² + (1/(2μ))||Λ^{-α}f||₂².

Consequently, on every finite preterminal interval [0,T],

    ||u(t)||₂² + μ∫_0^t ||Λ^αu||₂² ds
      ≤ ||u0||₂² + μ^{-1}∫_0^t ||Λ^{-α}f||₂² ds,   t≤T.

In particular, the sum of the supremum of the first term and the full
dissipation integral is at most twice the right side at T.

Write a(t)=||Λ^{5/4}u(t)||₂. Since α≥5/4 and the smallest nonzero
frequency is fixed, a²≤C(L,α)||Λ^αu||₂². Thus ∫_0^T a(t)² dt is paid
by the preceding energy estimate, uniformly as T approaches any fixed
finite hypothetical maximal time. This is the integrable coefficient
needed for Grönwall.

## Higher energy at the endpoint α=5/4

Fix an integer m≥3, put X_m=||u||_{H^m}² and
D_m=||Λ^αu||_{H^m}. Differentiate the equation through order m and pair
with the same derivatives of u. The term in which all derivatives fall
on the transported velocity cancels by incompressibility. Every remaining
term at order 1≤j≤m is a componentwise integral of the form

    ∫ (∇^j u)(∇^{j1}u)(∇^{j2}u),
    1≤j1,j2≤j,   j1+j2=j+1.

Here each ∇^j denotes an array of all derivatives of that order; the
finite contraction constants are absorbed into C_m.

The periodic Sobolev embeddings used are finite-exponent endpoints:

    ||∇u||_{12/5} ≤ C||Λ^{5/4}u||₂,
    ||∇^j u||_{12} ≤ C||Λ^{5/4}∇^j u||₂.

For j≥2 and 1≤b≤j, the Gagliardo–Nirenberg interpolation on ∇u gives

    ||∇^b u||_{p_b}
      ≤ C ||∇u||_{12/5}^{1-θ_b} ||∇^j u||₂^{θ_b},
    θ_b=(b-1)/(j-1),   1/p_b=(1-θ_b)·5/12+θ_b/2.

Because θ_{j1}+θ_{j2}=1, Hölder gives

    ||(∇^{j1}u)(∇^{j2}u)||_{12/11}
      ≤ C ||∇u||_{12/5} ||∇^j u||₂.

For j=1 the same product estimate follows directly from Hölder, placing
one copy of ∇u in L^{12/5} and the other in L². Pairing with the L^{12}
factor and using α≥5/4 yields, after summing the finite derivative terms,

    |higher nonlinear work| ≤ C_m a(t) X_m^{1/2} D_m
      ≤ (μ/2)D_m² + C_m μ^{-1}a(t)² X_m.

This is the critical estimate. No L∞ embedding at the critical exponent
and no extra high-norm integrability hypothesis have been used.

## Grönwall and continuation

The force obeys |<f,u>_{H^m}|≤(X_m+||f||_{H^m}²)/2. Discarding ordinary
viscous dissipation and changing C_m if necessary gives

    X_m' + μD_m²
      ≤ (1+C_m μ^{-1}a(t)²) X_m + ||f||_{H^m}².

The coefficient is integrable by native energy, so

    X_m(t) ≤ (X_m(0)+∫_0^T ||f||_{H^m}² ds)
              exp(T+C_m μ^{-1}∫_0^T a(s)² ds),   t≤T.

The right side is finite and independent of the preterminal endpoint.
Local strong well-posedness and its H^m continuation criterion therefore
exclude a finite maximal time. Applying the estimate at every integer
m≥3 propagates spatial smoothness; the equation propagates time
regularity for smooth forcing. Standard difference-energy uniqueness is
applied within this strong class. This is the classical hyperviscous
result, at its existing scope.

## Why the previous interpolation did not close

Absorbing a term of size
||Λ^{s+1}u||₂ ||Λ^s u||₂² by interpolation between H^s and H^{s+α}
produces a power of X_s with exponent (3α-1)/(2α-1)>1. That scalar
inequality alone does not give global bounds by Grönwall, even at α=5/4.
It has been replaced by the integrable-coefficient estimate above.

The whole-space critical energy method and its higher-derivative product
estimates are developed in [Tao, Analysis & PDE 2 (2009), Section 2,
pp. 363–366](https://msp.org/apde/2009/2-3/apde-v2-n3-p04-p.pdf).
The periodic proof above spells out the fixed-torus and forcing adapters;
it does not import the paper's stronger logarithmic result.
