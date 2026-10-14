import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Linarith

/-! Independent scalar absorption control. This is not a Navier tail producer. -/
namespace Fluids.UnconditionalNavier

theorem orthogonal_absorption (low high total K q : ℝ)
    (hdecomp : total = low + high) (hnative : low ≤ K)
    (hgap : high ≤ q * total) (hq : q < 1) :
    total ≤ K / (1 - q) := by
  apply (le_div_iff (by linarith : 0 < 1 - q)).2
  nlinarith

end Fluids.UnconditionalNavier
