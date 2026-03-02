import numpy as np

from nswave.deterministic_npz import savez_deterministic, sha256_file


def test_savez_deterministic_is_byte_stable(tmp_path):
    arr_a = np.arange(12, dtype=np.float64).reshape(3, 4)
    arr_b = np.array([3, 1, 4, 1, 5], dtype=np.int64)

    p1 = tmp_path / "a.npz"
    p2 = tmp_path / "b.npz"

    savez_deterministic(p1, z_last=arr_b, a_first=arr_a)
    savez_deterministic(p2, z_last=arr_b, a_first=arr_a)

    assert p1.read_bytes() == p2.read_bytes()
    assert sha256_file(p1) == sha256_file(p2)
