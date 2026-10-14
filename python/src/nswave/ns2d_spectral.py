"""2D pseudo-spectral vorticity solver on T^2 with de-aliasing."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np

from nswave.operators import kgrid_2d, nyquist_mask, spectral_time_grid, validate_real_spectrum


Array = np.ndarray


def dealias_mask_2d(N: int, frac: float = 2 / 3) -> Array:
    """Return a 2D de-aliasing mask using the 2/3 rule."""
    kk = np.fft.fftfreq(N)
    kcut = frac * (N // 2)
    kx, ky = np.meshgrid(kk * N, kk * N, indexing="ij")
    # Strict cutoff: modes at N/3 alias back into the retained band.
    mask = (np.abs(kx) < kcut) & (np.abs(ky) < kcut)
    return mask.astype(float)


def make_deterministic_omega0_hat(
    N: int,
    *,
    L: float = 2 * np.pi,
    seed: int = 0,
    k_init: Optional[int] = None,
    dealias: bool = True,
) -> Array:
    """Deterministic real-valued initial condition in Fourier space."""
    rng = np.random.default_rng(seed)
    omega0 = rng.standard_normal((N, N))
    omega_hat = np.fft.fft2(omega0)
    omega_hat[0, 0] = 0.0

    if k_init is not None:
        kk = np.fft.fftfreq(N) * N
        kx_int, ky_int = np.meshgrid(kk, kk, indexing="ij")
        mask_init = (kx_int * kx_int + ky_int * ky_int) <= k_init**2
        omega_hat = omega_hat * mask_init

    if dealias:
        omega_hat = omega_hat * dealias_mask_2d(N)
    _ = L
    return omega_hat * nyquist_mask(N, 2)


def vorticity_to_velocity_hat(omega_hat: Array, kx: Array, ky: Array, k2: Array) -> Array:
    """Compute velocity Fourier coefficients from vorticity."""
    psi_hat = np.zeros_like(omega_hat)
    mask = k2 != 0
    psi_hat[mask] = omega_hat[mask] / k2[mask]

    u_hat = np.empty((2,) + omega_hat.shape, dtype=complex)
    u_hat[0] = 1j * ky * psi_hat
    u_hat[1] = -1j * kx * psi_hat
    return u_hat


def nonlinear_term_hat(
    omega_hat: Array,
    kx: Array,
    ky: Array,
    k2: Array,
    *,
    forcing_hat_fn: Optional[Callable[[float, Array, Array], Array]] = None,
    t: float = 0.0,
    mask: Optional[Array] = None,
) -> Array:
    """Compute -u·∇ω + forcing in Fourier space (with optional de-aliasing)."""
    u_hat = vorticity_to_velocity_hat(omega_hat, kx, ky, k2)
    u = np.fft.ifft2(u_hat[0]).real, np.fft.ifft2(u_hat[1]).real

    dwdx_hat = 1j * kx * omega_hat
    dwdy_hat = 1j * ky * omega_hat
    dwdx = np.fft.ifft2(dwdx_hat).real
    dwdy = np.fft.ifft2(dwdy_hat).real

    adv = u[0] * dwdx + u[1] * dwdy
    adv_hat = np.fft.fft2(adv)

    if mask is not None:
        adv_hat = adv_hat * mask

    if forcing_hat_fn is None:
        return -adv_hat

    f_hat = np.asarray(forcing_hat_fn(t, kx, ky))
    if f_hat.shape != omega_hat.shape:
        raise ValueError("forcing spectrum shape must match vorticity")
    validate_real_spectrum(f_hat, (0, 1))
    if abs(f_hat[0, 0]) > 1e-12 * max(float(np.max(np.abs(f_hat))), np.finfo(float).tiny):
        raise ValueError("periodic vorticity forcing must have zero mean")
    f_hat = f_hat.copy()
    f_hat[0, 0] = 0
    if mask is not None:
        f_hat = f_hat * mask
    return -adv_hat + f_hat


def step_etdrk4(
    omega_hat: Array,
    dt: float,
    E: Array,
    E2: Array,
    linop: Array,
    nonlinear_fn: Callable[[Array, float], Array],
    t: float,
) -> Array:
    """Integrating-factor RK4 step; historical name retained for callers.

    Apply classical RK4 to exp(-L*t)*omega. This is not an ETDRK4
    phi-function scheme. The unused linop argument is kept for compatibility.
    """

    a = nonlinear_fn(omega_hat, t)
    b = nonlinear_fn(E2 * (omega_hat + (dt / 2) * a), t + dt / 2)
    c = nonlinear_fn(E2 * omega_hat + (dt / 2) * b, t + dt / 2)
    d = nonlinear_fn(E * omega_hat + dt * E2 * c, t + dt)

    return E * omega_hat + (dt / 6.0) * (E * a + 2.0 * E2 * b + 2.0 * E2 * c + d)


@dataclass
class NS2DResult:
    t: Array
    energy: Array
    enstrophy: Array
    diss_nu: Array
    diss_mu: Array
    omega_hat_final: Array
    omega_hat_hist: Optional[Array]
    t_omega: Optional[Array]
    success: bool
    message: str


def simulate_ns2d(
    *,
    N: int,
    dt: float,
    t_max: float,
    nu: float,
    mu: float,
    alpha: float,
    omega0_hat: Optional[Array] = None,
    forcing_hat_fn: Optional[Callable[[float, Array, Array], Array]] = None,
    L: float = 2 * np.pi,
    dealias: bool = True,
    record_every: int = 1,
    record_omega_hat_every: Optional[int] = None,
) -> NS2DResult:
    """Simulate 2D vorticity dynamics with energy dissipation diagnostics."""
    if N <= 0:
        raise ValueError("N must be positive")
    if dt <= 0:
        raise ValueError("dt must be positive")
    if t_max <= 0:
        raise ValueError("t_max must be positive")
    if record_every <= 0:
        raise ValueError("record_every must be positive")
    if record_omega_hat_every is not None and record_omega_hat_every <= 0:
        raise ValueError("record_omega_hat_every must be positive when set")

    kx, ky, k2 = kgrid_2d(N, L=L)
    if not np.all(np.isfinite([nu, mu, alpha])):
        raise ValueError("viscosity coefficients and alpha must be finite")
    if alpha < 1.0:
        raise ValueError("alpha must be >= 1.0 for energy dissipation diagnostics")

    k2_alpha = np.where(k2 > 0, k2 ** alpha, 0.0)
    k2_alpha_energy = np.where(k2 > 0, k2 ** ((alpha - 1.0) / 2.0), 0.0)
    linop = -(nu * k2 + mu * k2_alpha)

    mask = dealias_mask_2d(N) if dealias else nyquist_mask(N, 2)

    if omega0_hat is None:
        omega0_hat = make_deterministic_omega0_hat(N, L=L, dealias=dealias)

    omega_hat = np.asarray(omega0_hat, dtype=complex).copy()
    if omega_hat.shape != (N, N):
        raise ValueError("omega0_hat must have shape (N,N)")
    validate_real_spectrum(omega_hat, (0, 1))
    if abs(omega_hat[0, 0]) > 1e-12 * max(float(np.max(np.abs(omega_hat))), np.finfo(float).tiny):
        raise ValueError("periodic vorticity must have zero mean")
    omega_hat[0, 0] = 0
    if mask is not None:
        omega_hat = omega_hat * mask

    times = spectral_time_grid(dt, t_max)
    nsteps = len(times) - 1

    nrec = nsteps // record_every + 1 + int(nsteps % record_every != 0)
    t_hist = np.zeros(nrec)
    E_hist = np.zeros(nrec)
    Z_hist = np.zeros(nrec)
    diss_nu_hist = np.zeros(nrec)
    diss_mu_hist = np.zeros(nrec)

    def record_state(idx: int, t: float, omega_hat_state: Array) -> None:
        u_hat = vorticity_to_velocity_hat(omega_hat_state, kx, ky, k2)
        u0 = np.fft.ifft2(u_hat[0]).real
        u1 = np.fft.ifft2(u_hat[1]).real
        omega = np.fft.ifft2(omega_hat_state).real

        energy = 0.5 * np.mean(u0 * u0 + u1 * u1) * L * L
        enstrophy = 0.5 * np.mean(omega * omega) * L * L

        int_omega2 = np.mean(omega * omega) * L * L
        diss_nu = nu * int_omega2

        zeta = np.fft.ifft2(k2_alpha_energy * omega_hat_state).real
        int_zeta2 = np.mean(zeta * zeta) * L * L
        diss_mu = mu * int_zeta2
        if not np.all(np.isfinite([energy, enstrophy, diss_nu, diss_mu])):
            raise FloatingPointError("nonfinite 2D energy diagnostics")

        t_hist[idx] = t
        E_hist[idx] = energy
        Z_hist[idx] = enstrophy
        diss_nu_hist[idx] = diss_nu
        diss_mu_hist[idx] = diss_mu

    def nonlinear_fn(omega_hat_state: Array, t: float) -> Array:
        return nonlinear_term_hat(
            omega_hat_state,
            kx,
            ky,
            k2,
            forcing_hat_fn=forcing_hat_fn,
            t=t,
            mask=mask,
        )

    if record_omega_hat_every is None:
        omega_hist = None
        t_omega = None
    else:
        nrec_omega = nsteps // record_omega_hat_every + 1 + int(nsteps % record_omega_hat_every != 0)
        omega_hist = np.empty((nrec_omega, N, N), dtype=complex)
        t_omega = np.empty(nrec_omega)
        omega_hist[0] = omega_hat
        t_omega[0] = 0.0
        omega_rec_idx = 1

    record_state(0, 0.0, omega_hat)
    rec_idx = 1

    t = 0.0
    for step in range(1, nsteps + 1):
        dt_step = times[step] - times[step - 1]
        e_step = np.exp(linop * dt_step)
        e_half = np.exp(linop * dt_step / 2)
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            omega_hat = step_etdrk4(omega_hat, dt_step, e_step, e_half, linop, nonlinear_fn, t)
        if not np.all(np.isfinite(omega_hat)):
            raise FloatingPointError(f"nonfinite 2D state at step {step}")
        if mask is not None:
            omega_hat = omega_hat * mask
        t = float(times[step])

        if record_omega_hat_every is not None and (step % record_omega_hat_every == 0 or step == nsteps):
            omega_hist[omega_rec_idx] = omega_hat
            t_omega[omega_rec_idx] = t
            omega_rec_idx += 1

        if step % record_every == 0 or step == nsteps:
            record_state(rec_idx, t, omega_hat)
            rec_idx += 1

    return NS2DResult(
        t=t_hist,
        energy=E_hist,
        enstrophy=Z_hist,
        diss_nu=diss_nu_hist,
        diss_mu=diss_mu_hist,
        omega_hat_final=omega_hat,
        omega_hat_hist=omega_hist,
        t_omega=t_omega,
        success=True,
        message="ok",
    )
