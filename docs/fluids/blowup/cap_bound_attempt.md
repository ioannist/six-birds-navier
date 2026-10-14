# Bounds for the actual SGS supplied power

Use the sharp projections and supplied-force sign fixed in
[the bridge definition](ns_bridge_object.md): v=G_j u, w=(I-G_j)u,
y_j=Δ_jG_j(NL(v)-NL(u)), and p_j=<Δ_j u,y_j>.

## Correct energy bound

The tensor driving y_j is -(v⊗w+w⊗v+w⊗w). For a sharp shell projector
in three dimensions, the finite Fourier sum and Cauchy–Schwarz give

    ||Δ_jP∇·F||₂ ≤ C 2^{5j/2} ||F||₁.

The factor counts O(2^{3j}) frequencies with one derivative of size O(2^j).
It is valid for this finite-output projector; no L¹ boundedness of the
full Leray projection is asserted. Hence

    ||y_j||₂ ≤ C 2^{5j/2}(2||v||₂||w||₂+||w||₂²)
             ≤ C' 2^{5j/2} ||u||₂²,
    |p_j| ≤ C' 2^{5j/2} ||u||₂² ||u_j||₂.

For unforced strong NS with ||u(t)||₂≤M=||u0||₂ this gives the independent
**absolute** throughput bound |p_j|≤C'2^{5j/2}M³. This is a bound on power
per time, not an ICAP coefficient multiplying ||u_j||₂². Dividing the
first power bound by ||u_j||₂² introduces 1/||u_j||₂, which can be arbitrarily
large. If the shell is zero, its power is zero, but its response need not be.

The resulting exponential absolute capacity is a valid loose CAP choice.
For constant work quantum its reciprocal comparison series converges;
that choice fails to certify DIV. It gives neither convergence of the
reciprocal of the minimal capacity nor a general impossibility theorem.

## Exact counterexample to the previous shell-only bound

On the 2π torus write k=(4,0,0), l=(-3,1,0), q=k+l=(1,1,0),
a=(0,1,0), b=(1,3,0), c=(1,-1,0), and

    h(x)=a cos(k·x)+b cos(l·x),
    u_ε(x)=h(x)-ε c sin(q·x),   ε>0.

All three modes are divergence-free: k·a=l·b=q·c=0. They are smooth and
zero mean. At j=0, G_0u_ε=u_{ε,0}=-εc sin(q·x); both k and l are high modes.
The low mode has zero self-advection. The q-mode of NL(h) is
(3/2)c sin(q·x): the sum-interaction vector is
-(1/2)((a·l)b+(b·k)a)=-(1/2)(1,7,0), whose Leray projection at q is (3/2)c.
All mixed low/high interactions miss q and -q. Therefore

    y_0(u_ε)=-(3/2)c sin(q·x),
    p_0=(3/2) ε Vol(T³),
    ||u_{ε,0}||₂²=ε² Vol(T³).

The ratio p_0/||u_{ε,0}||₂²=3/(2ε) is unbounded, while ||∇u_ε||∞ stays
bounded. Consequently neither ||y_j||₂≲||∇u||∞||u_j||₂ nor
|p_j|≲||∇u||∞||u_j||₂² holds uniformly for this full SGS port. A smooth
local NS path exists for each such datum, so the control is in the original
instantaneous domain. The reproducible trigonometric test is
`python/tests/test_spectral_mathematics.py`.

## Scope of the correction

The earlier shell-only bound and its claimed uniform Λ(j) energy fallback
are withdrawn. Standard commutator estimates for selected high–low terms
cannot replace the full high–high source. No estimate on the nonlinear
packaged history class or all cumulative windows is proved here. The
abstract discrete diagonal-gain ICAP no-go in Lean remains valid.
