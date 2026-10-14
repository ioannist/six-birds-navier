import Fluids

-- Re-run with `lake env lean ScaffoldAxiomAudit.lean` from formal/.
-- These commands inspect the transitive theorem dependencies, not just source text.
#print axioms Fluids.inner_rhs_self
#print axioms Fluids.inner_rhs_self_le_zero
#print axioms Fluids.deriv_nonneg_of_nonneg_right
#print axioms Fluids.viscosity_coeff_nonneg
#print axioms Fluids.fixed_mode_passivity_counterexample
#print axioms Fluids.inner_rhsDiff_self
#print axioms Fluids.inner_rhsDiff_self_le_inner_g
#print axioms Fluids.hasDerivAt_normSq
#print axioms Fluids.deriv_normSq_le_two_inner_g
#print axioms Fluids.deriv_normSq_le_two_inner_g_at_rhs
#print axioms Fluids.le_t0_add_sum_range_of_increments_ge
#print axioms Fluids.monotone_of_increments_ge_of_nonneg
#print axioms Fluids.no_finite_upper_bound_of_increments_ge_of_sum_unbounded
#print axioms Fluids.no_finite_upper_bound_dyadic
#print axioms Fluids.workcap_div_le_increment
#print axioms Fluids.increments_ge_of_workcap
#print axioms Fluids.no_finite_upper_bound_of_workcap_of_sum_unbounded
#print axioms Fluids.no_finite_upper_bound_dyadic_of_workcap
#print axioms Fluids.telescoping_pow
#print axioms Fluids.route_capacity_poly_bound
#print axioms Fluids.route_capacity_poly_bound_of_bounded_sum
#print axioms Fluids.route_capacity_poly_bound_of_tame_factor
#print axioms Fluids.route_capacity_poly_bound_of_tame_factor_of_bounded_sum
#print axioms Fluids.max_add_le_add_max
#print axioms Fluids.max_sum_le_sum_max
#print axioms Fluids.icap_sum_le_card_mul
#print axioms Fluids.mass_cell_le_card_mul_deloc_mul_coeffEnergy
#print axioms Fluids.synthesis_energy_eq_coeffEnergy_of_orthonormal
#print axioms Fluids.mass_cell_le_card_mul_deloc_mul_fieldEnergy_of_orthonormal
#print axioms Fluids.icap_implies_pointwise_pos_bound_of_basis_tests
#print axioms Fluids.icap_implies_pointwise_pos_bound_on_feasible_basis_tests
