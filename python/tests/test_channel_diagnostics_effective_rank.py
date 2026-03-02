import numpy as np

from nswave.channel_diagnostics import effective_rank_m95_pooled, effective_rank_m_energy


def test_effective_rank_matches_known_rank_with_tiny_noise():
    rng = np.random.default_rng(7)
    n_samples = 40
    n_features = 24
    rank_true = 3

    A = rng.standard_normal((n_samples, rank_true))
    B = rng.standard_normal((rank_true, n_features))
    X = A @ B
    X = X + 1e-8 * rng.standard_normal(X.shape)

    svals = np.linalg.svd(X, full_matrices=False, compute_uv=False)
    m95 = effective_rank_m_energy(svals, 0.95)
    m95_pooled = effective_rank_m95_pooled(X, 0.95)
    assert abs(m95 - rank_true) <= 1
    assert m95_pooled == m95
