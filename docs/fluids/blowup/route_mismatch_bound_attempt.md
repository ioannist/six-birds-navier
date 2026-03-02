# Route mismatch bound attempt (NS filters)

## A. Setup and objects (NS terms; concise)
- Domain: T^3, divergence-free u(t,x), Leray projector 𝒫.
- Nonlinearity: NL(u) := 𝒫 ∇·(u⊗u).
- Filters: R_n = spectral truncation to level n (nested cutoff), NL_0(u) := R_0 NL(u).
- General filter form: τ_G(u) := G(u⊗u) − (Gu)⊗(Gu), C_G(u) := 𝒫∇·τ_G(u).

## B. Exact route mismatch identity (the algebra)
- Define u_0 := G_0 u, u_1 := G_1 u.
- Identity for nested/idempotent filters (G_0G_1 = G_0):
  τ_{G_0}(u) = τ_{G_0}(u_1) + G_0 τ_{G_1}(u).
- Mismatch tensor:
  Δτ := τ_{G_0}(u) − τ_{G_0}(u_1) − G_0 τ_{G_1}(u).
- For G_0G_1 = G_0, Δτ ≡ 0; hence ΔC := 𝒫∇·(Δτ) ≡ 0 and the shell‑projected mismatch is zero.

## C. Specialize to the repo ladder (spectral truncation)
- Identify G_0 ↔ R_n, G_1 ↔ R_{n+1}, and R_n R_{n+1} = R_n (nested spectral truncation).
- Conclusion: for the nested ladder, route mismatch is exactly zero at the field level.
- Numerical consistency: `python/scripts/nsbu21_route_capacity_numbers.py`, `python/scripts/nsbu23_route_capacity_hypothesis_check.py`.

## D. Non‑nested filters (optional general form)
- General mismatch for arbitrary G_0, G_1:
  Δτ = (G_0 − G_0G_1)(u⊗u) + (G_0G_1u⊗G_0G_1u − G_0u⊗G_0u).
- Any bound for ΔC := 𝒫∇·Δτ typically depends on a sup‑norm channel such as G(t)=||∇u(t)||_∞ unless additional structure is assumed.

## E. Binary conclusion block
- **Bound line:** For nested spectral truncations R_n, Δτ ≡ 0 (hence ΔC ≡ 0 and route mismatch vanishes).
- **Obstruction line:** For non‑nested filters, mismatch bounds typically depend on G(t)=||∇u(t)||_∞ (sup‑norm/vorticity channel).
