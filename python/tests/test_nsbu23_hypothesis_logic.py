import numpy as np

from nswave.route_capacity_hypothesis import fit_smallest_p, classify_summability, assess_mismatch_summability


def test_fit_smallest_p_good_case():
    steps = np.arange(6)
    r_steps = ((steps + 2) / (steps + 1)) ** 2
    p_best, max_violation = fit_smallest_p(r_steps, margin_rel=1e-6, p_max=5)
    assert p_best is not None
    assert p_best <= 2
    assert max_violation <= 1e-6


def test_fit_smallest_p_bad_case():
    r_steps = np.full(6, 1.8)
    p_best, _ = fit_smallest_p(r_steps, margin_rel=1e-3, p_max=3)
    assert p_best is None


def test_fit_smallest_p_offset():
    steps = np.arange(6)
    r_steps = ((steps + 3) / (steps + 2)) ** 2
    p_best, max_violation = fit_smallest_p(r_steps, margin_rel=1e-6, p_max=5, step_offset=1)
    assert p_best is not None
    assert p_best <= 2
    assert max_violation <= 1e-6


def test_classify_summability_pass():
    n = 40
    e = 1.0 / (np.arange(1, n + 1) ** 2)
    res = classify_summability(e, e, min_steps=8)
    assert res["status"] == "PASS"


def test_classify_summability_fail():
    e = np.ones(40)
    res = classify_summability(e, e, min_steps=8)
    assert res["status"] == "FAIL"


def test_classify_summability_inconclusive():
    n = 40
    e = 1.0 / (np.arange(1, n + 1))
    res = classify_summability(e, e, min_steps=8)
    assert res["status"] == "INCONCLUSIVE"


def test_assess_mismatch_trivial():
    e = np.zeros(5)
    e_tilde = np.zeros(5)
    lambda_direct = np.full(5, 1e-3)
    res = assess_mismatch_summability(e, e_tilde, lambda_direct)
    assert res["status"].startswith("PASS")
    assert res["q_est_raw"] is None
    assert res["tail_ratio_median_raw"] is None


def test_assess_mismatch_ignores_slack():
    e = np.zeros(5)
    e_tilde = np.zeros(5)
    lambda_direct = np.full(5, 1e-3)
    res = assess_mismatch_summability(e, e_tilde, lambda_direct)
    assert res["status"].startswith("PASS")
