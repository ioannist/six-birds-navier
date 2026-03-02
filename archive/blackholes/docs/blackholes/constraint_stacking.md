# Constraint stacking (toy)

- Z family: Debye-sum K=2, Z(omega) = a0 + a1/(1 - i omega/b1) + a2/(1 - i omega/b2).
- Fixed parameters: a0=0.2, a2=0.3, b2=5.0. Scan a1 in [0,1.2], b1 in [0.2,5.0].
- Barrier: Regge-Wheeler (M=1, ell=2), sampled on tortoise grid with peak shifted to x=0.
- Constraints:
  - Ringdown proxy: D_ring = max_{omega in [1,3]} |t(omega)^2 * R0(omega)| <= 0.05
  - Love-like proxy: L_proxy = |Re Z(0) - 1| <= 0.2
  - Low-omega absorption: A_low = 1 - |R0(omega_low)|^2 >= 0.20 (omega_low ~ 0.3)
  - No gain: G_low = max_{omega in [0.2,1.0]} |R0(omega) * r(omega)| <= 0.98
- Figure shows the number of constraints satisfied across the (a1, b1) grid, with feasible region outlined.
