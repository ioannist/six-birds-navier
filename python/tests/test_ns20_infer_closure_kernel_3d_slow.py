import numpy as np
import pytest

from nswave.coarse_grain import restrict_hat_3d, shell_index_from_k2
from nswave.ns3d_spectral import (
    make_deterministic_u0_hat_3d,
    nonlinear_term_hat_3d,
    simulate_ns3d,
)
from nswave.operators import dealias_mask_3d, kgrid_3d


@pytest.mark.slow
def test_ns20_closure_kernel_3d_slow() -> None:
    L = 2 * np.pi
    N_hi = 24
    N_co = 12
    dt = 5e-4
    t_max = 0.02
    nu = 1e-3
    mu = 0.0
    alpha = 2.0
    dealias = True
    record_u_hat_every = 10
    t_min = 0.005
    k_init = 8

    u0_hat = make_deterministic_u0_hat_3d(
        N_hi, L=L, seed=0, k_init=k_init, dealias=dealias
    )

    res_hi = simulate_ns3d(
        N=N_hi,
        dt=dt,
        t_max=t_max,
        nu=nu,
        mu=mu,
        alpha=alpha,
        u0_hat=u0_hat,
        forcing_hat_fn=None,
        L=L,
        dealias=dealias,
        record_every=10,
        record_u_hat_every=record_u_hat_every,
    )

    assert res_hi.u_hat_hist is not None
    assert res_hi.t_hat is not None

    kx_hi, ky_hi, kz_hi, k2_hi = kgrid_3d(N_hi, L=L)
    mask_hi = dealias_mask_3d(N_hi) if dealias else None

    kx, ky, kz, k2 = kgrid_3d(N_co, L=L)
    mask_co = dealias_mask_3d(N_co) if dealias else None

    shell = shell_index_from_k2(k2)
    s_max = int(shell.max())
    num = np.zeros(s_max + 1, dtype=float)
    den = np.zeros(s_max + 1, dtype=float)

    used = False
    for t, u_hi in zip(res_hi.t_hat, res_hi.u_hat_hist):
        if t < t_min:
            continue
        if mask_hi is not None:
            u_hi = u_hi * mask_hi

        u_co = restrict_hat_3d(u_hi, N_co)
        if mask_co is not None:
            u_co = u_co * mask_co

        lin_hi = -(nu * k2_hi)
        nl_hi = nonlinear_term_hat_3d(
            u_hi, kx_hi, ky_hi, kz_hi, k2_hi, mask=mask_hi, t=t, forcing_hat_fn=None
        )
        rhs_hi = lin_hi * u_hi + nl_hi
        rhs_hi_co = restrict_hat_3d(rhs_hi, N_co)

        lin_co = -(nu * k2)
        nl_co = nonlinear_term_hat_3d(
            u_co, kx, ky, kz, k2, mask=mask_co, t=t, forcing_hat_fn=None
        )
        rhs_no = lin_co * u_co + nl_co

        tau = rhs_hi_co - rhs_no
        if mask_co is not None:
            tau = tau * mask_co

        dot = np.real(
            tau[0] * np.conj(u_co[0])
            + tau[1] * np.conj(u_co[1])
            + tau[2] * np.conj(u_co[2])
        )
        uu = np.abs(u_co[0]) ** 2 + np.abs(u_co[1]) ** 2 + np.abs(u_co[2]) ** 2

        num += np.bincount(shell.ravel(), weights=-dot.ravel(), minlength=s_max + 1)
        den += np.bincount(shell.ravel(), weights=uu.ravel(), minlength=s_max + 1)
        used = True
        break

    assert used

    ell = num / np.maximum(den, 1e-30)
    dissip_total = float(np.sum(num[1:]))

    shells = np.arange(len(ell))
    low_mask = (shells >= 1) & (shells <= 3) & (den > 0)
    assert np.count_nonzero(low_mask) >= 3

    s_low = shells[low_mask]
    ell_low = ell[low_mask]
    nu_eff = float(np.sum(s_low**2 * ell_low) / np.sum(s_low**4))

    assert nu_eff > 1e-10
    assert dissip_total > -1e-8
