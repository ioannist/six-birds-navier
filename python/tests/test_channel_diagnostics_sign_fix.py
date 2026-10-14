import numpy as np

from nswave.channel_diagnostics import pca_channels


def test_pca_channels_sign_fix_is_deterministic_under_global_flip():
    rng = np.random.default_rng(11)
    X = rng.standard_normal((32, 20))
    # Break near-degeneracies to stabilize ordering in this test.
    X[:, 0] *= 3.0
    X[:, 1] *= 2.0
    X[:, 2] *= 1.5

    _, V1 = pca_channels(X, k=5)
    _, V2 = pca_channels(-X, k=5)

    assert V1.shape == V2.shape
    assert np.allclose(V1, V2, atol=1e-10, rtol=0.0)
    for row in V1:
        idx = int(np.argmax(np.abs(row)))
        assert row[idx] >= 0.0
