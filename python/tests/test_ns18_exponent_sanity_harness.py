import numpy as np

from nswave.operators import kgrid_3d, project_div_free_3d


def make_shell_field_u_hat(N: int, K: int, *, L: float, seed: int) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    u = rng.standard_normal((3, N, N, N))
    u_hat = np.fft.fftn(u, axes=(1, 2, 3))
    u_hat[:, 0, 0, 0] = 0.0

    kx, ky, kz, k2 = kgrid_3d(N, L=L)
    kk = np.fft.fftfreq(N) * N
    kx_i, ky_i, kz_i = np.meshgrid(kk, kk, kk, indexing="ij")
    kmag = np.sqrt(kx_i * kx_i + ky_i * ky_i + kz_i * kz_i)
    shell_mask = np.abs(kmag - K) <= 0.5

    u_hat = u_hat * shell_mask
    u_hat = project_div_free_3d(u_hat, kx, ky, kz)
    return u_hat, k2


def norm_lambda(u_hat: np.ndarray, k2: np.ndarray, r: float, L: float) -> float:
    mult = np.zeros_like(k2)
    mask = k2 > 0
    mult[mask] = k2[mask] ** (r / 2.0)
    u0 = np.fft.ifftn(mult * u_hat[0]).real
    u1 = np.fft.ifftn(mult * u_hat[1]).real
    u2 = np.fft.ifftn(mult * u_hat[2]).real
    val = np.mean(u0 * u0 + u1 * u1 + u2 * u2) * (L**3)
    return float(np.sqrt(max(val, 0.0)))


def slope_two_points(K1: float, K2: float, R1: float, R2: float) -> float:
    return (np.log(R2) - np.log(R1)) / (np.log(K2) - np.log(K1))


def test_ns18_interpolation_slope_separates_good_bad() -> None:
    N = 24
    L = 2 * np.pi
    s = 3.0
    alpha = 1.25
    theta_good = 1.0 / alpha
    theta_bad = 0.4
    K1, K2 = 4, 8

    u1, k2_1 = make_shell_field_u_hat(N, K1, L=L, seed=0)
    u2, k2_2 = make_shell_field_u_hat(N, K2, L=L, seed=0)

    A1 = norm_lambda(u1, k2_1, s, L)
    A1_1 = norm_lambda(u1, k2_1, s + 1.0, L)
    B1 = norm_lambda(u1, k2_1, s + alpha, L)
    R1_good = A1_1 / max((A1 ** (1.0 - theta_good)) * (B1**theta_good), 1e-30)
    R1_bad = A1_1 / max((A1 ** (1.0 - theta_bad)) * (B1**theta_bad), 1e-30)

    A2 = norm_lambda(u2, k2_2, s, L)
    A2_1 = norm_lambda(u2, k2_2, s + 1.0, L)
    B2 = norm_lambda(u2, k2_2, s + alpha, L)
    R2_good = A2_1 / max((A2 ** (1.0 - theta_good)) * (B2**theta_good), 1e-30)
    R2_bad = A2_1 / max((A2 ** (1.0 - theta_bad)) * (B2**theta_bad), 1e-30)

    slope_good = slope_two_points(K1, K2, R1_good, R2_good)
    slope_bad = slope_two_points(K1, K2, R1_bad, R2_bad)

    assert abs(slope_good) <= 0.20
    assert abs(slope_bad) >= 0.20
