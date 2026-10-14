import Mathlib.Analysis.Calculus.Deriv.Basic
import Mathlib.Analysis.Calculus.Deriv.Comp
import Mathlib.Order.Filter.Basic
import Mathlib.Tactic

namespace Fluids

open Filter
open scoped Topology

/-- If `f` is nonnegative to the right and differentiable at 0, then its derivative is nonnegative. -/
theorem deriv_nonneg_of_nonneg_right
    {f : ℝ → ℝ} {ν : ℝ}
    (h0 : f 0 = 0)
    (hder : HasDerivAt f ν 0)
    (hpos : ∀ᶠ x in 𝓝[>] (0:ℝ), 0 ≤ f x) :
    0 ≤ ν := by
  have hlim : Tendsto (fun t => t⁻¹ * f t) (𝓝[>] (0:ℝ)) (𝓝 ν) := by
    simpa [h0] using (hder.tendsto_slope_zero_right :
      Tendsto (fun t => t⁻¹ • (f (0 + t) - f 0)) (𝓝[>] 0) (𝓝 ν))
  have hposx : ∀ᶠ x in 𝓝[>] (0:ℝ), 0 < x := by
    simpa using (eventually_mem_nhdsWithin :
      ∀ᶠ x in 𝓝[Set.Ioi (0:ℝ)] (0:ℝ), x ∈ Set.Ioi (0:ℝ))
  have hpos' : ∀ᶠ x in 𝓝[>] (0:ℝ), 0 ≤ x⁻¹ * f x := by
    filter_upwards [hpos, hposx] with x hx hxpos
    have hxnonneg : 0 ≤ x := le_of_lt hxpos
    have hxinv : 0 ≤ x⁻¹ := inv_nonneg.mpr hxnonneg
    exact mul_nonneg hxinv hx
  by_contra hneg
  have hupper : ∀ᶠ x in 𝓝[>] (0:ℝ), (fun t => t⁻¹ * f t) x < 0 :=
    (tendsto_order.1 hlim).2 0 (by linarith)
  have hfalse : ∀ᶠ _ in 𝓝[>] (0:ℝ), False :=
    (hpos'.and hupper).mono (by
      intro _ hx
      exact (not_lt_of_ge hx.1) hx.2)
  have hbot : (𝓝[>] (0:ℝ)) = ⊥ :=
    (eventually_false_iff_eq_bot).1 hfalse
  exact (by
    have hne : (𝓝[>] (0:ℝ)).NeBot := by infer_instance
    exact hne.ne hbot)

/-- Viscosity coefficient is nonnegative under right-sided passivity. -/
theorem viscosity_coeff_nonneg
    {ell : ℝ → ℝ} {ν : ℝ}
    (h0 : ell 0 = 0)
    (hder : HasDerivAt ell ν 0)
    (hpass : ∀ᶠ r in 𝓝[>] (0:ℝ), 0 ≤ ell r) :
    0 ≤ ν :=
  deriv_nonneg_of_nonneg_right h0 hder hpass

/-- A fixed box's mode passivity (r >= 1) cannot replace right-neighborhood
passivity: this smooth symbol has a strictly negative derivative at zero. -/
theorem fixed_mode_passivity_counterexample :
    ∃ ell : ℝ → ℝ, ell 0 = 0 ∧ (∀ r : ℝ, 1 ≤ r → 0 ≤ ell r) ∧
      HasDerivAt ell (-(1 / 2 : ℝ)) 0 := by
  refine ⟨fun r => r ^ 2 - r / 2, ?_, ?_, ?_⟩
  · norm_num
  · intro r hr
    change 0 ≤ r ^ 2 - r / 2
    have hprod : 0 ≤ r * (r - 1 / 2) :=
      mul_nonneg (by linarith) (by linarith)
    nlinarith
  · simpa using ((hasDerivAt_id (0 : ℝ)).pow 2).sub
      ((hasDerivAt_id (0 : ℝ)).div_const 2)

end Fluids
