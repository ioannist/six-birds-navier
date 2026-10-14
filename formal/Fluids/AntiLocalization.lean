import Mathlib.Algebra.BigOperators.Group.Finset
import Mathlib.Algebra.Order.Chebyshev
import Mathlib.Data.Real.Basic
import Mathlib.Tactic

open scoped BigOperators

namespace Fluids

/-- Discrete anti-localization from finite delocalized channels. -/
theorem mass_cell_le_card_mul_deloc_mul_coeffEnergy
    {ι ρ : Type} [Fintype ι] [DecidableEq ι] [Fintype ρ] [DecidableEq ρ]
    (B : Finset ι)
    (ψ : ρ → ι → ℝ) (a : ρ → ℝ)
    (κ μB : ℝ)
    (hκ : 0 ≤ κ) (hμB : 0 ≤ μB)
    (hdeloc : ∀ r, (∑ i in B, (ψ r i) ^ 2) ≤ κ * μB) :
    (∑ i in B, ((∑ r, a r * ψ r i) ^ 2))
      ≤ (Fintype.card ρ : ℝ) * κ * μB * (∑ r, (a r) ^ 2) := by
  have _hκμB_nonneg : 0 ≤ κ * μB := mul_nonneg hκ hμB
  have hpoint :
      (∑ i in B, ((∑ r, a r * ψ r i) ^ 2))
        ≤ ∑ i in B, ((Fintype.card ρ : ℝ) * ∑ r, (a r * ψ r i) ^ 2) := by
    refine Finset.sum_le_sum ?_
    intro i _
    have hsq :
        (∑ r : ρ, a r * ψ r i) ^ 2
          ≤ (Fintype.card ρ : ℝ) * ∑ r : ρ, (a r * ψ r i) ^ 2 := by
      simpa [Finset.card_univ] using
        (sq_sum_le_card_mul_sum_sq (s := (Finset.univ : Finset ρ))
          (f := fun r : ρ => a r * ψ r i))
    exact hsq
  have hsplit :
      (∑ i in B, ∑ r : ρ, (a r * ψ r i) ^ 2)
        ≤ (κ * μB) * ∑ r : ρ, (a r) ^ 2 := by
    calc
      (∑ i in B, ∑ r : ρ, (a r * ψ r i) ^ 2)
          = ∑ r : ρ, ∑ i in B, (a r * ψ r i) ^ 2 := by
              rw [Finset.sum_comm]
      _ = ∑ r : ρ, ((a r) ^ 2 * ∑ i in B, (ψ r i) ^ 2) := by
            refine Finset.sum_congr rfl ?_
            intro r _
            calc
              (∑ i in B, (a r * ψ r i) ^ 2)
                  = ∑ i in B, ((a r) ^ 2 * (ψ r i) ^ 2) := by
                      refine Finset.sum_congr rfl ?_
                      intro i _
                      ring
              _ = (a r) ^ 2 * ∑ i in B, (ψ r i) ^ 2 := by
                    rw [Finset.mul_sum]
      _ ≤ ∑ r : ρ, ((a r) ^ 2 * (κ * μB)) := by
            refine Finset.sum_le_sum ?_
            intro r _
            have ha_nonneg : 0 ≤ (a r) ^ 2 := sq_nonneg (a r)
            exact mul_le_mul_of_nonneg_left (hdeloc r) ha_nonneg
      _ = (κ * μB) * ∑ r : ρ, (a r) ^ 2 := by
            rw [Finset.mul_sum]
            refine Finset.sum_congr rfl ?_
            intro r _
            ring
  calc
    (∑ i in B, ((∑ r, a r * ψ r i) ^ 2))
        ≤ ∑ i in B, ((Fintype.card ρ : ℝ) * ∑ r, (a r * ψ r i) ^ 2) := hpoint
    _ = (Fintype.card ρ : ℝ) * (∑ i in B, ∑ r, (a r * ψ r i) ^ 2) := by
          rw [Finset.mul_sum]
    _ ≤ (Fintype.card ρ : ℝ) * ((κ * μB) * ∑ r, (a r) ^ 2) := by
          have hcard_nonneg : 0 ≤ (Fintype.card ρ : ℝ) := by exact_mod_cast (Nat.zero_le _)
          exact mul_le_mul_of_nonneg_left hsplit hcard_nonneg
    _ = (Fintype.card ρ : ℝ) * κ * μB * (∑ r, (a r) ^ 2) := by ring

/-- Discrete Parseval identity from the channel Gram matrix, without assuming
the desired energy identity as a separate conclusion-carrying premise. -/
theorem synthesis_energy_eq_coeffEnergy_of_orthonormal
    {ι ρ : Type} [Fintype ι] [Fintype ρ] [DecidableEq ρ]
    (ψ : ρ → ι → ℝ) (a : ρ → ℝ)
    (horth : ∀ r s, (∑ i, ψ r i * ψ s i) = if r = s then 1 else 0) :
    (∑ i, (∑ r, a r * ψ r i) ^ 2) = ∑ r, (a r) ^ 2 := by
  classical
  calc
    (∑ i, (∑ r, a r * ψ r i) ^ 2)
        = ∑ i, ∑ r, ∑ s, (a r * a s) * (ψ r i * ψ s i) := by
            apply Finset.sum_congr rfl
            intro i _
            simp only [pow_two, Finset.sum_mul, Finset.mul_sum]
            apply Finset.sum_congr rfl
            intro r _
            apply Finset.sum_congr rfl
            intro s _
            ring
    _ = ∑ r, ∑ i, ∑ s, (a r * a s) * (ψ r i * ψ s i) := by
          rw [Finset.sum_comm]
    _ = ∑ r, ∑ s, ∑ i, (a r * a s) * (ψ r i * ψ s i) := by
          apply Finset.sum_congr rfl
          intro r _
          rw [Finset.sum_comm]
    _ = ∑ r, ∑ s, (a r * a s) * (∑ i, ψ r i * ψ s i) := by
          simp_rw [Finset.mul_sum]
    _ = ∑ r, (a r) ^ 2 := by
          simp [horth, pow_two]

/-- The field-energy corollary stated in the paper, with orthonormality explicit.
The delocalization hypothesis must hold for each cell B to obtain a uniform bound. -/
theorem mass_cell_le_card_mul_deloc_mul_fieldEnergy_of_orthonormal
    {ι ρ : Type} [Fintype ι] [DecidableEq ι] [Fintype ρ] [DecidableEq ρ]
    (B : Finset ι) (ψ : ρ → ι → ℝ) (a : ρ → ℝ) (κ μB : ℝ)
    (hκ : 0 ≤ κ) (hμB : 0 ≤ μB)
    (hdeloc : ∀ r, (∑ i in B, (ψ r i) ^ 2) ≤ κ * μB)
    (horth : ∀ r s, (∑ i, ψ r i * ψ s i) = if r = s then 1 else 0) :
    (∑ i in B, (∑ r, a r * ψ r i) ^ 2)
      ≤ (Fintype.card ρ : ℝ) * κ * μB * (∑ i, (∑ r, a r * ψ r i) ^ 2) := by
  rw [synthesis_energy_eq_coeffEnergy_of_orthonormal ψ a horth]
  exact mass_cell_le_card_mul_deloc_mul_coeffEnergy B ψ a κ μB hκ hμB hdeloc

end Fluids
