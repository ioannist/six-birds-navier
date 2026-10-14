import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Linarith

/-! Scalar spectral sign law only. The matrix/PDE/endpoint adapters are
written separately; this declaration does not pay an original NS source. -/
namespace Fluids.UnconditionalNavier

theorem middle_eigenvalue_production (a b c : ℝ)
    (hab : a ≤ b) (hbc : b ≤ c) (htrace : a + b + c = 0) :
    -4 * a * b * c ≤ 2 * max b 0 * (a ^ 2 + b ^ 2 + c ^ 2) := by
  have ha : a ≤ 0 := by linarith
  have hc : 0 ≤ c := by linarith
  have hac : 0 ≤ -a * c := mul_nonneg (neg_nonneg.mpr ha) hc
  have hprod : -a * c ≤ (a ^ 2 + c ^ 2) / 2 := by
    nlinarith [sq_nonneg (a + c)]
  by_cases hb : 0 ≤ b
  · rw [max_eq_left hb]
    have hm := mul_le_mul_of_nonneg_right hprod hb
    have hs := mul_nonneg hb (sq_nonneg b)
    nlinarith
  · have hb0 : b ≤ 0 := le_of_not_ge hb
    rw [max_eq_right hb0]
    have hm := mul_nonneg hac (neg_nonneg.mpr hb0)
    nlinarith

end Fluids.UnconditionalNavier
