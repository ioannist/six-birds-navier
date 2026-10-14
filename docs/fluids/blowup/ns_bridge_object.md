# Navier Stokes packaged bridge definition

The scaffold uses a nonlinear causal bridge on resolved histories, conditional
on a unique strong high-mode IVP. The sign and the projection convention below
make its readout the actual closure forcing in the resolved equation.
There is no claim of global existence for every candidate history.

## Fixed projections and spaces

Work on a fixed torus with zero-mean solenoidal H^s fields, s>5/2, ν>0,
and no external force. Let Δ_j be the **sharp orthogonal Fourier projector**
to S_j={k:2^j≤|k|<2^{j+1}}. Let G_j project to |k|<2^{j+1}, and H_j=I-G_j.
For zero-mean fields these shells decompose u. On a general torus the
lowest block must also include the nonzero frequencies below the first
shell; the diagnostics use side length 2π, whose first frequency is 1.

Thus Δ_j²=Δ_j, G_j²=G_j, Δ_jG_j=Δ_j, G_jH_j=0 and all these maps commute
with Leray P and spatial derivatives. U_j=Ran(Δ_j)∩H and
U_{≤j}=Ran(G_j)∩H. These are sharp spectral shells, not smooth overlapping
Littlewood–Paley multipliers. A smooth-filter formulation would require a
different retained/environment adapter and cannot use these identities.

Set NL(u)=P∇·(u⊗u). A strong actual solution obeys

    ∂t u = νΔu - NL(u).

## Sign of the actual resolved forcing

For v=G_j u and w=H_j u, define

    C_{≤j}(u) = G_j NL(v) - G_j NL(u),
    y_j(u) = Δ_j C_{≤j}(u).

Then the exact equations are

    ∂t v = νΔv - G_j NL(v) + C_{≤j}(u),
    ∂t u_j = νΔu_j - Δ_j NL(v) + y_j(u).

Therefore <u_j,y_j> is the closure's supplied power in the shell energy
identity. Reversing this sign exchanges supplied and extracted power and
changes positive work; it is not an innocuous notation choice. The same
sign is used in capacity_definition_ns.md and route_objects_ns.md.

## High-mode IVP and the history map

Fix u0∈H^s, T>0 and initial high state w0=H_j u0. Admit resolved histories
v∈C([0,T];U_{≤j}) for which the high-mode IVP has a unique strong solution
w∈C([0,T];H^s)∩L²([0,T];H^{s+1}), with H_j w=w:

    ∂t w = νΔw - H_j NL(v+w),   w(0)=w0.

If initial compatibility is part of the intended physical history class,
also impose v(0)=G_j u0. The high-mode range is invariant: applying G_j
leaves a homogeneous heat equation with zero initial data. Define

    Z_j[v](t) = Δ_j G_j (NL(v(t))-NL(v(t)+w(t))).

On a fixed torus this lies in U_j and is a continuous, hence L²-in-time,
shell readout. The input to Z_j is the **entire resolved history v**, not
just u_j; the power pairs its output with Δ_j v. Other resolved shells
can change the response while leaving Δ_j v unchanged.

Strong uniqueness also holds on every initial subinterval. Indeed, for
two high solutions with the same history, their difference z satisfies
an L² estimate controlled by the bounded strong gradients of v,w1,w2;
the transported z term cancels, and high projection drops from pairing
with z. Grönwall gives z=0. Consequently two admitted histories agreeing
on [0,t] give identical high states and outputs there. Continuity lets
one replace agreement by a.e. agreement. This establishes causality on
the admitted domain, not IVP existence for arbitrary histories.

## Exact return on an actual strong path

For an actual NS solution on [0,T], v=G_j u and w=H_j u solve the displayed
IVP. Its uniqueness identifies the packaged w with the actual high state,
so Z_j[G_j u]=y_j(u). The old smooth-filter/fixed-point definition did not
justify this return; the sharp projection convention repairs it without
assuming any new PDE estimate.

Frequency-response fits are empirical finite-window LTI approximations.
They neither define this nonlinear bridge nor establish ICAP on its
full admissible history class.
