"""3D pseudo-spectral incompressible NS solver on T^3."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np

from nswave.operators import dealias_mask_3d, kgrid_3d, project_div_free_3d

Array = np.ndarray


def make_deterministic_u0_hat_3d(
    N: int,
    *,
    L: float = 2 * np.pi,
    seed: int = 0,
    k_init: Optional[int] = 8,
    dealias: bool = True,
) -> Array:
    """Deterministic real-valued initial velocity field in Fourier space."""
    rng = np.random.default_rng(seed)
    u0 = rng.standard_normal((3, N, N, N))
    u0_hat = np.empty_like(u0, dtype=complex)
    for i in range(3):
        u0_hat[i] = np.fft.fftn(u0[i])
    u0_hat[:, 0, 0, 0] = 0.0

    if k_init is not None:
        kk = np.fft.fftfreq(N) * N
        kx_int, ky_int, kz_int = np.meshgrid(kk, kk, kk, indexing="ij")
        mask_init = (kx_int * kx_int + ky_int * ky_int + kz_int * kz_int) <= k_init**2
        u0_hat = u0_hat * mask_init

    kx, ky, kz, _ = kgrid_3d(N, L=L)
    u0_hat = project_div_free_3d(u0_hat, kx, ky, kz)

    if dealias:
        mask = dealias_mask_3d(N)
        u0_hat = u0_hat * mask

    return u0_hat


def _nonlinear_hat(
    u_hat: Array,
    kx: Array,
    ky: Array,
    kz: Array,
    k2: Array,
    *,
    mask: Optional[Array],
    t: float,
    forcing_hat_fn: Optional[Callable[[float, Array, Array, Array], Array]],
) -> Array:
    u = [np.fft.ifftn(u_hat[i]).real for i in range(3)]

    du = []
    for i in range(3):
        du.append(
            (
                np.fft.ifftn(1j * kx * u_hat[i]).real,
                np.fft.ifftn(1j * ky * u_hat[i]).real,
                np.fft.ifftn(1j * kz * u_hat[i]).real,
            )
        )

    adv = []
    for i in range(3):
        adv_i = u[0] * du[i][0] + u[1] * du[i][1] + u[2] * du[i][2]
        adv.append(adv_i)

    adv_hat = np.empty_like(u_hat)
    for i in range(3):
        adv_hat[i] = np.fft.fftn(adv[i])

    if mask is not None:
        adv_hat = adv_hat * mask

    nl_hat = -project_div_free_3d(adv_hat, kx, ky, kz)

    if forcing_hat_fn is None:
        return nl_hat

    f_hat = forcing_hat_fn(t, kx, ky, kz)
    f_hat = project_div_free_3d(f_hat, kx, ky, kz)
    if mask is not None:
        f_hat = f_hat * mask

    return nl_hat + f_hat


def nonlinear_term_hat_3d(
    u_hat: Array,
    kx: Array,
    ky: Array,
    kz: Array,
    k2: Array,
    *,
    mask: Optional[Array] = None,
    t: float = 0.0,
    forcing_hat_fn: Optional[Callable[[float, Array, Array, Array], Array]] = None,
) -> Array:
    """Return -P[(u·∇)u] (+ forcing if provided) in Fourier space."""
    return _nonlinear_hat(
        u_hat,
        kx,
        ky,
        kz,
        k2,
        mask=mask,
        t=t,
        forcing_hat_fn=forcing_hat_fn,
    )


def _step_ifrk4(
    u_hat: Array,
    dt: float,
    E: Array,
    E2: Array,
    L_hat: Array,
    nonlinear_fn: Callable[[Array, float], Array],
    t: float,
) -> Array:
    n1 = nonlinear_fn(u_hat, t)
    a = E2 * (u_hat + 0.5 * dt * n1)
    n2 = nonlinear_fn(a, t + 0.5 * dt)
    b = E2 * (u_hat + 0.5 * dt * n2)
    n3 = nonlinear_fn(b, t + 0.5 * dt)
    c = E * (u_hat + dt * n3)
    n4 = nonlinear_fn(c, t + dt)

    return E * u_hat + (dt / 6.0) * (E * n1 + 2.0 * E2 * n2 + 2.0 * E2 * n3 + n4)


@dataclass
class NS3DResult:
    t: Array
    energy: Array
    diss_nu: Array
    diss_mu: Array
    div_rms: Array
    u_hat_final: Array
    u_hat_hist: Optional[Array]
    t_hat: Optional[Array]
    success: bool
    message: str


def simulate_ns3d(
    *,
    N: int,
    dt: float,
    t_max: float,
    nu: float,
    mu: float,
    alpha: float,
    u0_hat: Optional[Array] = None,
    forcing_hat_fn: Optional[Callable[[float, Array, Array, Array], Array]] = None,
    L: float = 2 * np.pi,
    dealias: bool = True,
    record_every: int = 1,
    record_u_hat_every: Optional[int] = None,
) -> NS3DResult:
    if N <= 0:
        raise ValueError("N must be positive")
    if dt <= 0:
        raise ValueError("dt must be positive")
    if t_max <= 0:
        raise ValueError("t_max must be positive")
    if record_every <= 0:
        raise ValueError("record_every must be positive")
    if record_u_hat_every is not None and record_u_hat_every <= 0:
        raise ValueError("record_u_hat_every must be positive when set")
    if alpha < 1.0:
        raise ValueError("alpha must be >= 1.0")

    kx, ky, kz, k2 = kgrid_3d(N, L=L)
    mask = dealias_mask_3d(N) if dealias else None

    if u0_hat is None:
        u0_hat = make_deterministic_u0_hat_3d(N, L=L, dealias=dealias)

    u_hat = project_div_free_3d(u0_hat, kx, ky, kz)
    if mask is not None:
        u_hat = u_hat * mask

    k2_alpha = np.where(k2 > 0, k2 ** alpha, 0.0)
    L_hat = -(nu * k2 + mu * k2_alpha)
    L_hat[0, 0, 0] = 0.0
    E = np.exp(L_hat * dt)
    E2 = np.exp(L_hat * (dt / 2.0))

    def nonlinear_fn(u_hat_state: Array, t: float) -> Array:
        return _nonlinear_hat(
            u_hat_state,
            kx,
            ky,
            kz,
            k2,
            mask=mask,
            t=t,
            forcing_hat_fn=forcing_hat_fn,
        )

    nsteps = int(round(t_max / dt))
    if nsteps < 1:
        raise ValueError("t_max must be at least one dt")

    nrec = nsteps // record_every + 1
    t_hist = np.zeros(nrec)
    E_hist = np.zeros(nrec)
    diss_nu_hist = np.zeros(nrec)
    diss_mu_hist = np.zeros(nrec)
    div_hist = np.zeros(nrec)

    def record_state(idx: int, t: float, u_hat_state: Array) -> None:
        u = [np.fft.ifftn(u_hat_state[i]).real for i in range(3)]
        energy = 0.5 * np.mean(u[0] ** 2 + u[1] ** 2 + u[2] ** 2) * L**3

        div = np.fft.ifftn(1j * (kx * u_hat_state[0] + ky * u_hat_state[1] + kz * u_hat_state[2])).real
        div_rms = np.sqrt(np.mean(div * div))

        grad_sq = 0.0
        for i in range(3):
            dux = np.fft.ifftn(1j * kx * u_hat_state[i]).real
            duy = np.fft.ifftn(1j * ky * u_hat_state[i]).real
            duz = np.fft.ifftn(1j * kz * u_hat_state[i]).real
            grad_sq += dux * dux + duy * duy + duz * duz
        diss_nu = nu * np.mean(grad_sq) * L**3

        w_hat = np.zeros_like(u_hat_state)
        for i in range(3):
            w_hat[i] = (k2 ** (alpha / 2.0)) * u_hat_state[i]
        w = [np.fft.ifftn(w_hat[i]).real for i in range(3)]
        diss_mu = mu * np.mean(w[0] ** 2 + w[1] ** 2 + w[2] ** 2) * L**3

        t_hist[idx] = t
        E_hist[idx] = energy
        diss_nu_hist[idx] = diss_nu
        diss_mu_hist[idx] = diss_mu
        div_hist[idx] = div_rms

    record_state(0, 0.0, u_hat)
    rec_idx = 1

    if record_u_hat_every is None:
        u_hat_hist = None
        t_hat = None
    else:
        nrec_hat = nsteps // record_u_hat_every + 1
        u_hat_hist = np.empty((nrec_hat, 3, N, N, N), dtype=complex)
        t_hat = np.empty(nrec_hat)
        u_hat_hist[0] = u_hat
        t_hat[0] = 0.0
        hat_idx = 1
    t = 0.0

    for step in range(1, nsteps + 1):
        u_hat = _step_ifrk4(u_hat, dt, E, E2, L_hat, nonlinear_fn, t)
        if mask is not None:
            u_hat = u_hat * mask
        t += dt

        if record_u_hat_every is not None and step % record_u_hat_every == 0:
            u_hat_hist[hat_idx] = u_hat
            t_hat[hat_idx] = t
            hat_idx += 1

        if step % record_every == 0:
            record_state(rec_idx, t, u_hat)
            rec_idx += 1

    return NS3DResult(
        t=t_hist,
        energy=E_hist,
        diss_nu=diss_nu_hist,
        diss_mu=diss_mu_hist,
        div_rms=div_hist,
        u_hat_final=u_hat,
        u_hat_hist=u_hat_hist,
        t_hat=t_hat,
        success=True,
        message="ok",
    )
