import Mathlib

open Complex

lemma one_add_ne_zero_of_re_nonneg (Z : Complex) (hZ : 0 ≤ Z.re) : (1 + Z) ≠ 0 := by
  intro h
  have h' : Z = -1 := by
    have h'' : Z + 1 = 0 := by
      simpa [add_comm] using h
    exact eq_neg_of_add_eq_zero_left h''
  have : (0 : Real) ≤ (-1 : Real) := by
    simpa [h'] using hZ
  linarith

/-- Cayley-transform inequality for the impedance mapping. -/
theorem abs_R_from_Z_le_one_of_re_nonneg (Z : Complex) (hZ : 0 ≤ Z.re) :
    Complex.abs ((1 - Z) / (1 + Z)) ≤ 1 := by
  have hne : (1 + Z) ≠ 0 := one_add_ne_zero_of_re_nonneg Z hZ
  have habs : Complex.abs ((1 - Z) / (1 + Z)) =
      Complex.abs (1 - Z) / Complex.abs (1 + Z) := by
    simp
  have hsq : Complex.abs (1 - Z) ^ 2 ≤ Complex.abs (1 + Z) ^ 2 := by
    have h :
        (1 - Z.re) * (1 - Z.re) + Z.im * Z.im ≤
        (1 + Z.re) * (1 + Z.re) + Z.im * Z.im := by
      nlinarith [hZ]
    simpa [Complex.sq_abs, Complex.normSq_apply] using h
  have hineq : Complex.abs (1 - Z) ≤ Complex.abs (1 + Z) := by
    have h' : |Complex.abs (1 - Z)| ≤ |Complex.abs (1 + Z)| := (sq_le_sq).1 hsq
    simpa using h'
  have hpos : 0 < Complex.abs (1 + Z) := by
    exact Complex.abs.pos hne
  have hdiv : Complex.abs (1 - Z) / Complex.abs (1 + Z) ≤ 1 := by
    exact (div_le_one hpos).2 hineq
  simpa [habs] using hdiv

theorem re_nonneg_of_abs_cayley_le_one (Z : Complex)
    (hden : (1 + Z) ≠ 0)
    (h : Complex.abs ((1 - Z) / (1 + Z)) ≤ 1) : 0 ≤ Z.re := by
  have habs : Complex.abs ((1 - Z) / (1 + Z)) =
      Complex.abs (1 - Z) / Complex.abs (1 + Z) := by
    simp
  have hpos : 0 < Complex.abs (1 + Z) := by
    exact Complex.abs.pos hden
  have hdiv : Complex.abs (1 - Z) / Complex.abs (1 + Z) ≤ 1 := by
    simpa [habs] using h
  have hineq : Complex.abs (1 - Z) ≤ Complex.abs (1 + Z) := by
    exact (div_le_one hpos).1 hdiv
  have hsq : Complex.abs (1 - Z) ^ 2 ≤ Complex.abs (1 + Z) ^ 2 := by
    have h' : |Complex.abs (1 - Z)| ≤ |Complex.abs (1 + Z)| := by
      simpa using hineq
    exact (sq_le_sq).2 h'
  have hsq' :
      (1 - Z.re) * (1 - Z.re) + Z.im * Z.im ≤
      (1 + Z.re) * (1 + Z.re) + Z.im * Z.im := by
    simpa [Complex.sq_abs, Complex.normSq_apply] using hsq
  nlinarith [hsq']

theorem abs_cayley_le_one_iff_re_nonneg (Z : Complex) (hden : (1 + Z) ≠ 0) :
    Complex.abs ((1 - Z) / (1 + Z)) ≤ 1 ↔ 0 ≤ Z.re := by
  constructor
  · intro h
    exact re_nonneg_of_abs_cayley_le_one Z hden h
  · intro hZ
    exact abs_R_from_Z_le_one_of_re_nonneg Z hZ

-- BH-T2 hook
theorem passive_imp_absR_le_one (Z : Complex) (hZ : 0 ≤ Z.re) :
    Complex.abs ((1 - Z) / (1 + Z)) ≤ 1 :=
  abs_R_from_Z_le_one_of_re_nonneg Z hZ
