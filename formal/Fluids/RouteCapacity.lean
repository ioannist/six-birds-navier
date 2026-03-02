import Mathlib.Algebra.BigOperators.Group.Finset
import Mathlib.Data.Real.Basic
import Mathlib.Tactic

open scoped BigOperators

namespace Fluids

lemma telescoping_pow (n p : ℕ) :
    ((n + 2 : ℝ) / (n + 1)) ^ p * (n + 1) ^ p = (n + 2) ^ p := by
  have hb : (n + 1 : ℝ) ≠ 0 := by
    exact_mod_cast Nat.succ_ne_zero n
  field_simp [div_pow, hb]

theorem route_capacity_poly_bound
    {C e : ℕ → ℝ} (p : ℕ)
    (hC0 : 0 ≤ C 0)
    (he : ∀ n, 0 ≤ e n)
    (hrec : ∀ n, C (n + 1) ≤ ((n + 2 : ℝ) / (n + 1)) ^ p * C n + e n) :
    ∀ n, C n ≤ (C 0 + ∑ k in Finset.range n, e k) * (n + 1) ^ p := by
  intro n
  induction n with
  | zero =>
      simp [hC0]
  | succ n ih =>
      let r : ℝ := ((n + 2 : ℝ) / (n + 1)) ^ p
      have hcoef_nonneg : 0 ≤ r := by
        have hnum : (0 : ℝ) ≤ (n + 2 : ℝ) := by
          nlinarith
        have hden : (0 : ℝ) < (n + 1 : ℝ) := by
          nlinarith
        have hbase : 0 ≤ (n + 2 : ℝ) / (n + 1) := by
          exact div_nonneg hnum (le_of_lt hden)
        simpa [r] using (pow_nonneg hbase p)
      have hmul_le :
          r * C n ≤ r * ((C 0 + ∑ k in Finset.range n, e k) * (n + 1) ^ p) := by
        exact mul_le_mul_of_nonneg_left ih hcoef_nonneg
      have hrec' :
          C (n + 1) ≤ r * ((C 0 + ∑ k in Finset.range n, e k) * (n + 1) ^ p) + e n := by
        have hrecn := hrec n
        linarith [hrecn, hmul_le]
      have hmul :
          r * ((C 0 + ∑ k in Finset.range n, e k) * (n + 1) ^ p)
            = (C 0 + ∑ k in Finset.range n, e k) * (n + 2) ^ p := by
        have htel : r * (n + 1) ^ p = (n + 2) ^ p := by
          simpa [r] using (telescoping_pow n p)
        calc
          r * ((C 0 + ∑ k in Finset.range n, e k) * (n + 1) ^ p)
              = (C 0 + ∑ k in Finset.range n, e k) * (r * (n + 1) ^ p) := by
                ac_rfl
          _ = (C 0 + ∑ k in Finset.range n, e k) * (n + 2) ^ p := by
                simp [htel]
      have hstep :
          C (n + 1) ≤ (C 0 + ∑ k in Finset.range n, e k) * (n + 2) ^ p + e n := by
        simpa [hmul] using hrec'
      have hpow : (1 : ℝ) ≤ (n + 2 : ℝ) ^ p := by
        have hbase : (1 : ℝ) ≤ (n + 2 : ℝ) := by
          nlinarith
        have hnonneg : (0 : ℝ) ≤ (1 : ℝ) := by
          nlinarith
        have hpow' : (1 : ℝ) ^ p ≤ (n + 2 : ℝ) ^ p := by
          exact pow_le_pow_left hnonneg hbase p
        simpa using hpow'
      have hsum : (C 0 + ∑ k in Finset.range n, e k) * (n + 2) ^ p + e n
          ≤ (C 0 + ∑ k in Finset.range (n + 1), e k) * (n + 2) ^ p := by
        have hsum' : ∑ k in Finset.range (n + 1), e k
            = ∑ k in Finset.range n, e k + e n := by
              simp [Finset.sum_range_succ, add_comm, add_left_comm, add_assoc]
        have hpow' : e n ≤ e n * (n + 2 : ℝ) ^ p := by
          have hen : 0 ≤ e n := he n
          nlinarith [hpow, hen]
        calc
          (C 0 + ∑ k in Finset.range n, e k) * (n + 2) ^ p + e n
              ≤ (C 0 + ∑ k in Finset.range n, e k) * (n + 2) ^ p + e n * (n + 2 : ℝ) ^ p := by
                nlinarith [hpow']
          _ = (C 0 + ∑ k in Finset.range n, e k + e n) * (n + 2) ^ p := by
                ring
          _ = (C 0 + ∑ k in Finset.range (n + 1), e k) * (n + 2) ^ p := by
                nlinarith [hsum']
      have hfinal' : C (n + 1) ≤ (C 0 + ∑ k in Finset.range (n + 1), e k) * (n + 2) ^ p :=
        le_trans hstep hsum
      have hcast : (n + 2 : ℝ) = (n + 1 + 1 : ℝ) := by
        nlinarith
      have hfinal : C (n + 1) ≤ (C 0 + ∑ k in Finset.range (n + 1), e k) * (n + 1 + 1 : ℝ) ^ p := by
        simpa [hcast] using hfinal'
      simpa using hfinal


theorem route_capacity_poly_bound_of_bounded_sum
    {C e : ℕ → ℝ} (p : ℕ) {E : ℝ}
    (hC0 : 0 ≤ C 0)
    (he : ∀ n, 0 ≤ e n)
    (hrec : ∀ n, C (n + 1) ≤ ((n + 2 : ℝ) / (n + 1)) ^ p * C n + e n)
    (hsum : ∀ n, (∑ k in Finset.range n, e k) ≤ E) :
    ∀ n, C n ≤ (C 0 + E) * (n + 1) ^ p := by
  intro n
  have hbase := route_capacity_poly_bound (C:=C) (e:=e) p hC0 he hrec n
  have hsum' : C 0 + ∑ k in Finset.range n, e k ≤ C 0 + E := by
    nlinarith [hsum n]
  have hpow : 0 ≤ (n + 1 : ℝ) ^ p := by
    exact pow_nonneg (by nlinarith : (0 : ℝ) ≤ (n + 1 : ℝ)) p
  have hmul := mul_le_mul_of_nonneg_right hsum' hpow
  exact le_trans hbase hmul

end Fluids
