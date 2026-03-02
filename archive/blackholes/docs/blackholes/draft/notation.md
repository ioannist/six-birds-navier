This file is canonical for the BH section's notation.
All BH scripts/figures assume these conventions.
Changing it requires updating both paper and code.

Canonical definitions
- Time convention: time dependence exp(-i omega t)
- Coordinate: y := x - x0
- Plane-wave decomposition near a free region:
  psi(y,omega) = A_in(omega) e^{-i omega y} + A_out(omega) e^{+i omega y}
  R(omega) := A_out / A_in
- Admittance and impedance at inner boundary x=x0:
  Y(omega) := psi'(x0,omega) / psi(x0,omega)
  Z(omega) := i Y(omega) / omega
  equivalently: psi'(x0,omega) = - i omega Z(omega) psi(x0,omega)
- Mapping identities:
  R(Z) = (1 - Z) / (1 + Z)
  Z(R) = (1 - R) / (1 + R)
  R(Y) = (i omega + Y) / (i omega - Y)
  Y(R) = - i omega (1 - R) / (1 + R)
  Y(Z) = - i omega Z
  Z(Y) = i Y / omega
- Flux/current convention used for passivity:
  j := Im(psi* psi')
  for psi = e^{-i omega y} + R e^{+i omega y}, j = omega(|R|^2 - 1)

Implementation anchors
- docs/blackholes/conventions.md
- python/src/bhwave/conventions.py
- python/scripts/bh02_convention_sanity.py
