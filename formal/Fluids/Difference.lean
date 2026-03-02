import Mathlib.Analysis.InnerProductSpace.Basic
import Mathlib.Analysis.Calculus.Deriv.Basic
import Mathlib.Analysis.InnerProductSpace.Calculus
import Mathlib.Tactic

namespace Fluids

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]

/-- Difference RHS with linear dissipation and forcing. -/
def rhsDiff (A Aα : E →ₗ[ℝ] E) (nu mu : ℝ) (w g : E) : E :=
  -(nu • (A w)) - (mu • (Aα w)) + g

/-- Inner-product identity for the difference RHS. -/
theorem inner_rhsDiff_self (A Aα : E →ₗ[ℝ] E) (nu mu : ℝ) (w g : E) :
    ⟪rhsDiff A Aα nu mu w g, w⟫_ℝ
      = -nu * ⟪A w, w⟫_ℝ - mu * ⟪Aα w, w⟫_ℝ + ⟪g, w⟫_ℝ := by
  simp [rhsDiff, sub_eq_add_neg, inner_add_left, inner_add_right, inner_smul_left]

/-- With PSD operators and nonnegative coefficients, ⟪rhsDiff w, w⟫ ≤ ⟪g, w⟫. -/
theorem inner_rhsDiff_self_le_inner_g
    (A Aα : E →ₗ[ℝ] E) (nu mu : ℝ)
    (hA : ∀ w : E, 0 ≤ ⟪A w, w⟫_ℝ)
    (hAα : ∀ w : E, 0 ≤ ⟪Aα w, w⟫_ℝ)
    (hnu : 0 ≤ nu) (hmu : 0 ≤ mu) (w g : E) :
    ⟪rhsDiff A Aα nu mu w g, w⟫_ℝ ≤ ⟪g, w⟫_ℝ := by
  have hId := inner_rhsDiff_self (A:=A) (Aα:=Aα) (nu:=nu) (mu:=mu) (w:=w) (g:=g)
  have hA' := hA w
  have hAα' := hAα w
  nlinarith [hId, hA', hAα', hnu, hmu]

/-- Derivative of the norm-square along a differentiable curve. -/
theorem hasDerivAt_normSq
    {w : ℝ → E} {w' : E} {t : ℝ} (hw : HasDerivAt w w' t) :
    HasDerivAt (fun τ => ‖w τ‖^2) (2 * ⟪w', w t⟫_ℝ) t := by
  -- Use the norm-square derivative and symmetry of the real inner product.
  simpa [real_inner_comm] using hw.norm_sq

/-- Differential inequality for the difference equation derivative value. -/
theorem deriv_normSq_le_two_inner_g
    (A Aα : E →ₗ[ℝ] E) (nu mu : ℝ)
    (hA : ∀ w : E, 0 ≤ ⟪A w, w⟫_ℝ)
    (hAα : ∀ w : E, 0 ≤ ⟪Aα w, w⟫_ℝ)
    (hnu : 0 ≤ nu) (hmu : 0 ≤ mu)
    {w g : ℝ → E} {w' : E} {t : ℝ}
    (_hw : HasDerivAt w w' t)
    (hw_rhs : w' = rhsDiff A Aα nu mu (w t) (g t)) :
    (2 * ⟪w', w t⟫_ℝ) ≤ 2 * ⟪g t, w t⟫_ℝ := by
  have h_inner := inner_rhsDiff_self_le_inner_g (A:=A) (Aα:=Aα) (nu:=nu) (mu:=mu)
      hA hAα hnu hmu (w t) (g t)
  have h_inner' : ⟪w', w t⟫_ℝ ≤ ⟪g t, w t⟫_ℝ := by
    simpa [hw_rhs] using h_inner
  nlinarith [h_inner']

end Fluids
