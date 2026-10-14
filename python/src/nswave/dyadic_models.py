"""Toy dyadic shell ODE models (diagnostic, not a proof object)."""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from nswave.operators import spectral_time_grid

Array = np.ndarray


def kp_rhs(
    u: Array,
    *,
    nu: float,
    mu: float,
    alpha: float,
    lam: float = 2.0,
    forcing: Array | None = None,
) -> Array:
    """Katz–Pavlović-style dyadic cascade RHS."""
    J = u.size
    f = np.zeros_like(u)
    if forcing is not None:
        f = forcing

    du = np.zeros_like(u)
    for j in range(J):
        u_jm1 = u[j - 1] if j - 1 >= 0 else 0.0
        u_j = u[j]
        u_jp1 = u[j + 1] if j + 1 < J else 0.0
        du[j] = (
            (lam ** j) * (u_jm1**2)
            - (lam ** (j + 1)) * u_j * u_jp1
            - nu * (lam ** (2 * j)) * u_j
            - mu * (lam ** (2 * alpha * j)) * u_j
            + f[j]
        )
    return du


def _rk4_step(u: Array, dt: float, rhs_fn) -> Array:
    k1 = rhs_fn(u)
    k2 = rhs_fn(u + 0.5 * dt * k1)
    k3 = rhs_fn(u + 0.5 * dt * k2)
    k4 = rhs_fn(u + dt * k3)
    return u + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)


def simulate_kp(
    *,
    J: int,
    dt: float,
    t_max: float,
    nu: float,
    mu: float,
    alpha: float,
    u0: Array,
    lam: float = 2.0,
    forcing: Array | None = None,
    record_every: int = 1,
) -> dict:
    """Simulate the KP dyadic model with RK4."""
    if not isinstance(J, (int, np.integer)) or J <= 0:
        raise ValueError("J must be positive")
    u0 = np.asarray(u0, dtype=float)
    if u0.shape != (J,) or not np.all(np.isfinite(u0)):
        raise ValueError("u0 must have length J")
    if not np.all(np.isfinite([nu, mu, alpha, lam])) or lam <= 0:
        raise ValueError("model coefficients must be finite and lam positive")
    # Signed nu/mu remain allowed for the deliberately active controls.
    if forcing is not None:
        forcing = np.asarray(forcing, dtype=float)
        if forcing.shape != (J,) or not np.all(np.isfinite(forcing)):
            raise ValueError("forcing must be a finite vector of length J")
    if not isinstance(record_every, (int, np.integer)) or record_every <= 0:
        raise ValueError("record_every must be a positive integer")
    times = spectral_time_grid(dt, t_max)
    nsteps = len(times) - 1
    nrec = nsteps // record_every + 1 + int(nsteps % record_every != 0)
    t_hist = np.zeros(nrec)
    u_hist = np.zeros((nrec, J))
    E_hist = np.zeros(nrec)
    P_hist = np.zeros(nrec)

    u = u0.copy()

    def rhs_fn(x: Array) -> Array:
        return kp_rhs(x, nu=nu, mu=mu, alpha=alpha, lam=lam, forcing=forcing)

    rec_idx = 0
    for step in range(nsteps + 1):
        if step % record_every == 0 or step == nsteps:
            t_hist[rec_idx] = times[step]
            u_hist[rec_idx] = u
            E_hist[rec_idx] = 0.5 * float(np.sum(u * u))
            P_hist[rec_idx] = mu * float(np.sum((lam ** (2 * alpha * np.arange(J))) * (u * u)))
            if not np.all(np.isfinite([E_hist[rec_idx], P_hist[rec_idx]])):
                return {
                    "t": t_hist[:rec_idx], "u_hist": u_hist[:rec_idx],
                    "E_hist": E_hist[:rec_idx], "sgs_power_hist": P_hist[:rec_idx],
                    "success": False, "message": "Nonfinite energy/power diagnostic",
                }
            rec_idx += 1
        if step < nsteps:
            u = _rk4_step(u, times[step + 1] - times[step], rhs_fn)
            if not np.all(np.isfinite(u)):
                return {
                    "t": t_hist[:rec_idx],
                    "u_hist": u_hist[:rec_idx],
                    "E_hist": E_hist[:rec_idx],
                    "sgs_power_hist": P_hist[:rec_idx],
                    "success": False,
                    "message": "NaN encountered",
                }

    return {
        "t": t_hist,
        "u_hist": u_hist,
        "E_hist": E_hist,
        "sgs_power_hist": P_hist,
        "success": True,
        "message": "OK",
    }
