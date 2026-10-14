import Mathlib.Algebra.BigOperators.Group.Finset
import Mathlib.Data.Real.Basic
import Mathlib.Tactic

open scoped BigOperators

namespace Fluids

lemma max_add_le_add_max (x y : ℝ) : max (x + y) 0 ≤ max x 0 + max y 0 := by
  by_cases hx : 0 ≤ x
  · have hx' : max x 0 = x := by simp [max_eq_left hx]
    by_cases hy : 0 ≤ y
    · have hy' : max y 0 = y := by simp [max_eq_left hy]
      have hxy : 0 ≤ x + y := add_nonneg hx hy
      simp [hx', hy', max_eq_left hxy, add_assoc, add_comm, add_left_comm]
    · have hy' : max y 0 = 0 := by simp [max_eq_right (le_of_not_ge hy)]
      have hxy_le : x + y ≤ x := by linarith
      have hxy_max_le : max (x + y) 0 ≤ x := by
        exact (max_le_iff).2 ⟨hxy_le, hx⟩
      simpa [hx', hy'] using hxy_max_le
  · have hx' : max x 0 = 0 := by simp [max_eq_right (le_of_not_ge hx)]
    by_cases hy : 0 ≤ y
    · have hy' : max y 0 = y := by simp [max_eq_left hy]
      have hxy_le : x + y ≤ y := by linarith
      have hxy_max_le : max (x + y) 0 ≤ y := by
        exact (max_le_iff).2 ⟨hxy_le, hy⟩
      simpa [hx', hy'] using hxy_max_le
    · have hy' : max y 0 = 0 := by simp [max_eq_right (le_of_not_ge hy)]
      have hxy_nonpos : x + y ≤ 0 := by linarith
      have hxy' : max (x + y) 0 = 0 := by simp [max_eq_right hxy_nonpos]
      simp [hx', hy', hxy']


theorem max_sum_le_sum_max {ι : Type} (s : Finset ι) (a : ι → ℝ) :
    max (∑ i in s, a i) 0 ≤ ∑ i in s, max (a i) 0 := by
  classical
  refine Finset.induction_on s ?base ?step
  · simp
  · intro i s hi hs
    have hxy := max_add_le_add_max (a i) (∑ j in s, a j)
    calc
      max (∑ j in insert i s, a j) 0
          = max (a i + ∑ j in s, a j) 0 := by
              simp [Finset.sum_insert, hi, add_comm, add_left_comm, add_assoc]
      _ ≤ max (a i) 0 + max (∑ j in s, a j) 0 := hxy
      _ ≤ max (a i) 0 + ∑ j in s, max (a j) 0 := by
            exact add_le_add_left hs _
      _ = ∑ j in insert i s, max (a j) 0 := by
            simp [Finset.sum_insert, hi, add_comm, add_left_comm, add_assoc]


-- ICAP aggregation: if each atom has constant Λ0, the sum has constant m * Λ0.
theorem icap_sum_le_card_mul {ι : Type} (s : Finset ι)
    (Λ0 E : ℝ) (w : ι → ℝ)
    (h : ∀ i ∈ s, w i ≤ Λ0 * E) :
    ∑ i in s, w i ≤ (s.card : ℝ) * Λ0 * E := by
  classical
  revert h
  refine Finset.induction_on s ?base ?step
  · intro _
    simp
  · intro i s hi hs h
    have hi_bound : w i ≤ Λ0 * E := h i (by simp [hi])
    have hs_bound : ∑ j in s, w j ≤ (s.card : ℝ) * Λ0 * E := by
      apply hs
      intro j hj
      exact h j (by simp [hj, hi])
    have hsum : ∑ j in insert i s, w j = w i + ∑ j in s, w j := by
      simp [Finset.sum_insert, hi, add_comm, add_left_comm, add_assoc]
    have hcard : (insert i s).card = s.card + 1 := by
      exact Finset.card_insert_of_not_mem hi
    calc
      ∑ j in insert i s, w j
          = w i + ∑ j in s, w j := hsum
      _ ≤ Λ0 * E + (s.card : ℝ) * Λ0 * E := by
            linarith [hi_bound, hs_bound]
      _ = ((s.card : ℝ) + 1) * Λ0 * E := by ring
      _ = ((insert i s).card : ℝ) * Λ0 * E := by
            simp [hcard, Nat.cast_add, Nat.cast_one]

end Fluids
