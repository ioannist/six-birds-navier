import Mathlib.Analysis.InnerProductSpace.Basic
import Mathlib.LinearAlgebra.Basic
import Mathlib.Tactic

namespace Fluids

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]

/-- Abstract Navier-Stokes-style RHS with quadratic cancellation and linear dissipation. -/
def rhs (B : E → E → E) (A Aα : E →ₗ[ℝ] E) (nu mu : ℝ) (u : E) : E :=
  -(B u u) - (nu • (A u)) - (mu • (Aα u))

/-- Nonlinear cancellation leaves only dissipative contributions in ⟪rhs u, u⟫. -/
theorem inner_rhs_self (B : E → E → E) (A Aα : E →ₗ[ℝ] E) (nu mu : ℝ)
    (hB : ∀ u : E, ⟪B u u, u⟫_ℝ = 0) (u : E) :
    ⟪rhs B A Aα nu mu u, u⟫_ℝ
      = -nu * ⟪A u, u⟫_ℝ - mu * ⟪Aα u, u⟫_ℝ := by
  have hB' : ⟪B u u, u⟫_ℝ = 0 := hB u
  -- Expand inner products using linearity and the cancellation hypothesis.
  simp [
    rhs,
    sub_eq_add_neg,
    hB',
    inner_add_left,
    inner_add_right,
    inner_smul_left,
    inner_smul_right,
    mul_comm,
    mul_left_comm,
    mul_assoc,
  ]

/-- With PSD operators and nonnegative coefficients, ⟪rhs u, u⟫ is nonpositive. -/
theorem inner_rhs_self_le_zero (B : E → E → E) (A Aα : E →ₗ[ℝ] E) (nu mu : ℝ)
    (hB : ∀ u : E, ⟪B u u, u⟫_ℝ = 0)
    (hA : ∀ u : E, 0 ≤ ⟪A u, u⟫_ℝ)
    (hAα : ∀ u : E, 0 ≤ ⟪Aα u, u⟫_ℝ)
    (hnu : 0 ≤ nu) (hmu : 0 ≤ mu) (u : E) :
    ⟪rhs B A Aα nu mu u, u⟫_ℝ ≤ 0 := by
  have hId := inner_rhs_self (B:=B) (A:=A) (Aα:=Aα) (nu:=nu) (mu:=mu) hB u
  have hA' := hA u
  have hAα' := hAα u
  -- Combine the identity with nonnegativity of the dissipative forms.
  nlinarith [hId, hA', hAα', hnu, hmu]

end Fluids
