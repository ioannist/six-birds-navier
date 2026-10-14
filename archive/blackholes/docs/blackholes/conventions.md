# BH-02 conventions (internal)

## Time convention
- Physical time dependence uses exp(-i omega t).

## Coordinates and boundary location
- Exterior coordinate x increases outward.
- Inner boundary is at x = x0; define y := x - x0 so the boundary is y = 0 and the exterior is y >= 0.

## Near-boundary wave decomposition
- In a free-wave zone (V ~= 0), write:
  psi(y, omega) = A_in(omega) e^{-i omega y} + A_out(omega) e^{+i omega y}.
- Define the reflection coefficient as R(omega) := A_out / A_in.
- If one instead uses global x in the exponentials, R gains a trivial phase:
  R_x(omega) = e^{-2 i omega x0} R(omega). We do not use this as the primary definition.

## Boundary condition and response variables
- Normalized impedance Z(omega) (dimensionless) is defined by:
  psi'(0, omega) = - i omega Z(omega) psi(0, omega).
- Dimensional ratio Y(omega) := psi'(0, omega) / psi(0, omega), so Y(omega) = - i omega Z(omega).

## Mapping between R, Z, and Y
- R(omega) = (1 - Z(omega)) / (1 + Z(omega))
- Z(omega) = (1 - R(omega)) / (1 + R(omega))
- R(omega) = (i omega + Y(omega)) / (i omega - Y(omega))
- Y(omega) = - i omega (1 - R(omega)) / (1 + R(omega))
- Derivation: psi(0) = A_in + A_out, psi'(0) = -i omega A_in + i omega A_out.
  With R = A_out / A_in, Y = psi'/psi gives Y = -i omega (1 - R)/(1 + R).

## Flux and passivity
- Use the Schr/Wronskian current j := Im(psi* psi') (flux in +x).
- For real omega > 0 and the decomposition above:
  j = omega (|A_out|^2 - |A_in|^2).
- A passive boundary (no energy injection into the exterior) means j(0) <= 0 for an incident wave,
  which is equivalent to |R(omega)| <= 1.
- Algebraic equivalence with Z:
  |R|^2 = |1 - Z|^2 / |1 + Z|^2
       = ((1 - Re Z)^2 + (Im Z)^2) / ((1 + Re Z)^2 + (Im Z)^2).
  Thus |R| <= 1 <=> Re Z >= 0 (for real omega and no superradiance).

## Sanity check table
- Perfect absorber (purely ingoing): Z = 1 => R = 0; Y = -i omega.
- Neumann reflector: Z = 0 => R = 1; Y = 0.
- Dirichlet/hard wall limit: Z -> infinity => R -> -1.
- Lossless/reactive: Re Z = 0 => |R| = 1.

## Relation to Step-3 doc notation
Some Step-3 notes write psi' = Z psi with a different sign convention. In this repo:
Y = psi'/psi and Z_norm = -Y/(i omega). Overall phase conventions can flip R by a
global phase; we fix the definitions above as the standard.
