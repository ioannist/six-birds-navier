import numpy as np

from nswave.operators import (
    hypervisc_multiplier,
    kgrid_2d,
    kgrid_3d,
    laplacian_multiplier,
    project_div_free_2d,
    project_div_free_3d,
)


def test_project_div_free_2d_kills_divergence():
    N = 16
    rng = np.random.default_rng(0)
    u_hat = rng.standard_normal((2, N, N)) + 1j * rng.standard_normal((2, N, N))

    kx, ky, k2 = kgrid_2d(N)
    u_proj = project_div_free_2d(u_hat, kx, ky)

    div_hat = 1j * (kx * u_proj[0] + ky * u_proj[1])
    mask = k2 != 0
    scale = np.max(np.abs(u_hat))
    if scale == 0:
        scale = 1.0
    assert np.max(np.abs(div_hat[mask])) < 1e-10 * scale


def test_project_div_free_2d_idempotent():
    N = 16
    rng = np.random.default_rng(0)
    u_hat = rng.standard_normal((2, N, N)) + 1j * rng.standard_normal((2, N, N))

    kx, ky, _ = kgrid_2d(N)
    p1 = project_div_free_2d(u_hat, kx, ky)
    p2 = project_div_free_2d(p1, kx, ky)

    denom = np.linalg.norm(p1.ravel())
    if denom == 0:
        denom = 1.0
    err = np.linalg.norm((p2 - p1).ravel()) / denom
    assert err < 1e-12


def test_project_div_free_3d_kills_divergence_and_idempotent():
    N = 8
    rng = np.random.default_rng(0)
    u_hat = rng.standard_normal((3, N, N, N)) + 1j * rng.standard_normal((3, N, N, N))

    kx, ky, kz, k2 = kgrid_3d(N)
    u_proj = project_div_free_3d(u_hat, kx, ky, kz)

    div_hat = 1j * (kx * u_proj[0] + ky * u_proj[1] + kz * u_proj[2])
    mask = k2 != 0
    scale = np.max(np.abs(u_hat))
    if scale == 0:
        scale = 1.0
    assert np.max(np.abs(div_hat[mask])) < 1e-10 * scale

    p2 = project_div_free_3d(u_proj, kx, ky, kz)
    denom = np.linalg.norm(u_proj.ravel())
    if denom == 0:
        denom = 1.0
    err = np.linalg.norm((p2 - u_proj).ravel()) / denom
    assert err < 1e-12


def test_multipliers():
    N = 16
    _, _, k2 = kgrid_2d(N)

    hv = hypervisc_multiplier(k2, alpha=1.0)
    lap = laplacian_multiplier(k2)

    assert np.allclose(hv, k2)
    assert np.allclose(lap, -k2)
