import numpy as np

from nswave.channel_diagnostics import (
    channel_span_snapshot_bounds, effective_rank_m95_pooled, effective_rank_m_energy,
)


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


def test_rank_is_amplitude_invariant_and_total_energy_rank_stays_in_range():
    X = np.diag([3.0, 2.0, 1.0])
    for amplitude in [1.0, 1e-30, 1e-150]:
        assert effective_rank_m95_pooled(amplitude * X, 0.95) == 3
    assert effective_rank_m_energy([1.0, 1.0, 1.0], 1.0) == 3


def test_small_concentrated_residual_must_be_charged_in_channel_return():
    # Unit constant channel, and an orthogonal residual concentrated in cell 0.
    channel = np.ones(4)/2
    residual = np.array([np.sqrt(3)/2, -1/(2*np.sqrt(3)), -1/(2*np.sqrt(3)), -1/(2*np.sqrt(3))])
    snapshot = np.sqrt(0.95)*channel + np.sqrt(0.05)*residual
    local_fraction = snapshot[0]**2 / np.sum(snapshot**2)
    assert local_fraction > 0.25  # The span-only bound fails on this snapshot.
    for scale in [1.0, 1e-150]:
        delta, bounds = channel_span_snapshot_bounds((scale*snapshot)[None], channel[None], 0.25)
        assert np.isclose(delta[0], 0.05)
        assert bounds[0] >= local_fraction
        assert np.isclose(bounds[0], (np.sqrt(0.25*0.95)+np.sqrt(0.05))**2)
