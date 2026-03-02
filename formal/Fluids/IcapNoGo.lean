import Mathlib.Algebra.BigOperators.Group.Finset
import Mathlib.Data.Real.Basic
import Mathlib.Tactic

open scoped BigOperators

namespace Fluids

theorem icap_implies_pointwise_pos_bound_of_basis_tests
    {ι : Type} [Fintype ι] [DecidableEq ι]
    (g : ι → ℝ) (Λ : ℝ)
    (hicap : ∀ u : ι → ℝ,
      max (∑ i, g i * (u i) ^ 2) 0 ≤ Λ * (∑ i, (u i) ^ 2)) :
    ∀ i, max (g i) 0 ≤ Λ := by
  intro i
  have htest := hicap (fun j => if j = i then (1 : ℝ) else 0)
  simpa using htest

end Fluids
