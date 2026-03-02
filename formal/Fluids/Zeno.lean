import Mathlib.Algebra.BigOperators.Group.Finset
import Mathlib.Data.Real.Basic
import Mathlib.Order.Monotone.Basic
import Mathlib.Tactic

open scoped BigOperators

namespace Fluids


theorem le_t0_add_sum_range_of_increments_ge
  {t Φ : ℕ → ℝ}
  (hinc : ∀ n, Φ n ≤ t (n + 1) - t n) :
  ∀ n, t 0 + ∑ k in Finset.range n, Φ k ≤ t n := by
  intro n
  induction n with
  | zero =>
      simp
  | succ n ih =>
      have hsum : t 0 + (∑ k in Finset.range n, Φ k) + Φ n ≤ t n + Φ n := by
        linarith [ih]
      have hstep : t n + Φ n ≤ t (n + 1) := by
        linarith [hinc n]
      calc
        t 0 + ∑ k in Finset.range (n + 1), Φ k
            = t 0 + (∑ k in Finset.range n, Φ k) + Φ n := by
                simp [Finset.sum_range_succ, add_assoc, add_left_comm, add_comm]
        _ ≤ t n + Φ n := hsum
        _ ≤ t (n + 1) := hstep


theorem monotone_of_increments_ge_of_nonneg
  {t Φ : ℕ → ℝ}
  (hΦ : ∀ n, 0 ≤ Φ n)
  (hinc : ∀ n, Φ n ≤ t (n + 1) - t n) :
  Monotone t := by
  have hstep : ∀ n, t n ≤ t (n + 1) := by
    intro n
    have hnonneg : 0 ≤ t (n + 1) - t n := by
      exact le_trans (hΦ n) (hinc n)
    exact sub_nonneg.mp hnonneg
  exact monotone_nat_of_le_succ hstep


theorem no_finite_upper_bound_of_increments_ge_of_sum_unbounded
  {t Φ : ℕ → ℝ}
  (hinc : ∀ n, Φ n ≤ t (n + 1) - t n)
  (hunb : ∀ M : ℝ, ∃ n : ℕ, M ≤ ∑ k in Finset.range n, Φ k) :
  ¬ (∃ T : ℝ, ∀ n : ℕ, t n ≤ T) := by
  intro hT
  rcases hT with ⟨T, hT⟩
  have hsum_le : ∀ n, ∑ k in Finset.range n, Φ k ≤ T - t 0 := by
    intro n
    have hlower := le_t0_add_sum_range_of_increments_ge hinc n
    have hupper := hT n
    linarith [hlower, hupper]
  rcases hunb (T - t 0 + 1) with ⟨n, hn⟩
  have := hsum_le n
  linarith


def Φdyad (Φ : ℕ → ℝ) (j : ℕ) : ℝ := Φ (2 ^ j)


theorem no_finite_upper_bound_dyadic
  {t Φ : ℕ → ℝ}
  (hinc : ∀ j, Φ (2 ^ j) ≤ t (j + 1) - t j)
  (hunb : ∀ M : ℝ, ∃ n : ℕ, M ≤ ∑ k in Finset.range n, Φ (2 ^ k)) :
  ¬ (∃ T : ℝ, ∀ j : ℕ, t j ≤ T) := by
  have hinc' : ∀ j, Φdyad Φ j ≤ t (j + 1) - t j := by
    intro j
    exact hinc j
  have hunb' : ∀ M : ℝ, ∃ n : ℕ, M ≤ ∑ k in Finset.range n, Φdyad Φ k := by
    intro M
    simpa [Φdyad] using hunb M
  exact no_finite_upper_bound_of_increments_ge_of_sum_unbounded hinc' hunb'


end Fluids
