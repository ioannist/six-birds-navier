import numpy as np

from nswave.dyadic import crossing_times, monotone_envelope


def test_latency_extraction_monotone() -> None:
    t = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    j_series = np.array([0, 1, 0, 2, 2])
    j_mon = monotone_envelope(j_series)
    j_vals, t_reach = crossing_times(t, j_mon)

    assert np.array_equal(j_mon, np.array([0, 1, 1, 2, 2]))
    assert np.array_equal(j_vals, np.array([0, 1, 2]))
    assert np.allclose(t_reach, np.array([0.0, 1.0, 3.0]))
