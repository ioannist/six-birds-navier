import Fluids.Zeno
import Mathlib.Algebra.BigOperators.Group.Finset
import Mathlib.Data.Real.Basic
import Mathlib.Tactic

open scoped BigOperators

namespace Fluids

theorem workcap_div_le_increment
    {t w Cap : ℕ → ℝ} {n : ℕ}
    (hCap : 0 < Cap n)
    (h : w n ≤ Cap n * (t (n + 1) - t n)) :
    w n / Cap n ≤ t (n + 1) - t n := by
  have h' : w n ≤ (t (n + 1) - t n) * Cap n := by
    simpa [mul_comm] using h
  have hdiv := (div_le_iff hCap).2 h'
  simpa using hdiv


theorem increments_ge_of_workcap
    {t w Cap : ℕ → ℝ}
    (hCap : ∀ n, 0 < Cap n)
    (hworkcap : ∀ n, w n ≤ Cap n * (t (n + 1) - t n)) :
    ∀ n, w n / Cap n ≤ t (n + 1) - t n := by
  intro n
  exact workcap_div_le_increment (hCap n) (hworkcap n)


theorem no_finite_upper_bound_of_workcap_of_sum_unbounded
    {t w Cap : ℕ → ℝ}
    (hCap : ∀ n, 0 < Cap n)
    (hworkcap : ∀ n, w n ≤ Cap n * (t (n + 1) - t n))
    (hunb : ∀ M : ℝ, ∃ n : ℕ, M ≤ ∑ k in Finset.range n, (w k / Cap k)) :
    ¬ (∃ T : ℝ, ∀ n : ℕ, t n ≤ T) := by
  have hinc : ∀ n, (w n / Cap n) ≤ t (n + 1) - t n := by
    intro n
    exact workcap_div_le_increment (hCap n) (hworkcap n)
  exact no_finite_upper_bound_of_increments_ge_of_sum_unbounded hinc hunb


theorem no_finite_upper_bound_dyadic_of_workcap
    {t w Cap : ℕ → ℝ}
    (hCap : ∀ j, 0 < Cap (2 ^ j))
    (hworkcap : ∀ j, w (2 ^ j) ≤ Cap (2 ^ j) * (t (j + 1) - t j))
    (hunb : ∀ M : ℝ, ∃ n : ℕ, M ≤ ∑ k in Finset.range n, (w (2 ^ k) / Cap (2 ^ k))) :
    ¬ (∃ T : ℝ, ∀ j : ℕ, t j ≤ T) := by
  have hinc : ∀ j, (w (2 ^ j) / Cap (2 ^ j)) ≤ t (j + 1) - t j := by
    intro j
    have h' : w (2 ^ j) ≤ (t (j + 1) - t j) * Cap (2 ^ j) := by
      simpa [mul_comm] using hworkcap j
    have hdiv := (div_le_iff (hCap j)).2 h'
    simpa using hdiv
  exact no_finite_upper_bound_dyadic (t:=t) (Φ:=fun n => w n / Cap n) hinc hunb

end Fluids
