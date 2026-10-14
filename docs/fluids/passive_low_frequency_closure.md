# Passive low frequency closure at its justified scope

For a **linear, time-local scalar isotropic symbol** ℓ(|k|²) on
solenoidal velocities, suppose ℓ is C² near r=0, ℓ(0)=0, and ℓ(r)≥0
for all sufficiently small r>0. Then

    ℓ(r)=ηr+R(r),    η=ℓ'(0)≥0,    |R(r)|≤(M/2)r²,

where M bounds |ℓ''| on the interval. In a projected momentum equation
already containing νA, the leading coefficient is **ν+η**, and the
remaining closure term has a low-frequency symbol O(|k|⁴).

This is the existing low-frequency closure implication with its actual
premises and return made explicit. It does not construct a microscopic
closure or derive these premises from Six Birds role labels.

## Derivative and Taylor return

Since ℓ(0)=0, the difference quotient ℓ(r)/r is nonnegative for small
r>0. Differentiability at zero gives η≥0; this is precisely the
right-neighborhood passivity premise of `Fluids.viscosity_coeff_nonneg`.
The integral Taylor formula is

    R(r)=∫_0^r (r-s)ℓ''(s) ds,

so |R(r)|≤Mr²/2. For the linear scalar symbol, this identifies ηA as
the allowed order-two contribution. It is not a global high-frequency
bound or a regularity criterion.

η can vanish. The symbol ℓ(r)=r² is passive and smooth, vanishes on
constants, and has a fourth-order first nonzero term. Therefore a
**nonzero** leading viscosity is not forced without nondegeneracy. The
Laplacian is the first allowed term in this idealization, not always the
first nonzero one.

## The fixed torus gap

Passivity on a fixed periodic box constrains ℓ only at its realized
frequencies. It does not by itself constrain an interpolated scalar
symbol on a right neighborhood of zero. On the 2π torus, for example,

    ℓ(r)=r(r-1/2)

vanishes at zero and is strictly positive at every nonzero squared
Fourier frequency r=|k|²≥1, yet ℓ'(0)=-1/2. Thus even a smooth passive
scalar torus multiplier need not have a nonnegative derivative for an
arbitrarily chosen extension between zero and the first lattice point.

The derivative implication is valid when right-neighborhood positivity
is explicitly supplied, or is independently obtained from a constitutive
continuum-symbol law (for example by testing all small frequencies in a
continuum model or uniformly across a specified family of boxes). A
finite numerical low-k fit supplies neither such law nor an IR limit
on a fixed box. These are applicability assumptions of the existing
conditional result.

## Symmetry and memory assumptions

Translation invariance and self-adjointness alone do not imply a scalar
radial multiplier. The scalar isotropic constitutive idealization must
be explicit; full rotational covariance/parity arguments require their
own domain and representation assumptions. Rotations of the cubic
lattice are not all spatial rotations. This review does not add an
unproved symmetry classification theorem to the scaffold.

An autonomous lift of a memory equation is local in the **enlarged state**.
It does not make its velocity-only constitutive response time-local.
Consequently finite-memory autonomy is not a proof of the time-local
symbol premise above. The actual packaged SGS bridge retains its full
resolved history and initial high state.

Likewise, L² passivity or UV domination of a general nonlinear closure
does not automatically provide the differentiated coercivity and
commutator bounds used for a high Sobolev regularity theorem. The
hyperviscous dossier proves the chosen exact μA^α representative for
α≥5/4; it does not prove global regularity for every closure carrying
the abstract C0–C4 labels. No minimality or uniquely forced physical
completion follows from that example alone.
