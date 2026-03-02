"""Toy dyadic shell ODE models (diagnostic, not a proof object)."""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np

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
    if J <= 0:
        raise ValueError("J must be positive")
    if u0.size != J:
        raise ValueError("u0 must have length J")
    if dt <= 0 or t_max <= 0:
        raise ValueError("dt and t_max must be positive")

    nsteps = int(round(t_max / dt))
    if nsteps < 1:
        raise ValueError("t_max must be at least one dt")

    nrec = nsteps // record_every + 1
    t_hist = np.zeros(nrec)
    u_hist = np.zeros((nrec, J))
    E_hist = np.zeros(nrec)
    P_hist = np.zeros(nrec)

    u = u0.copy()

    def rhs_fn(x: Array) -> Array:
        return kp_rhs(x, nu=nu, mu=mu, alpha=alpha, lam=lam, forcing=forcing)

    rec_idx = 0
    t = 0.0
    for step in range(nsteps + 1):
        if step % record_every == 0:
            t_hist[rec_idx] = t
            u_hist[rec_idx] = u
            E_hist[rec_idx] = 0.5 * float(np.sum(u * u))
            P_hist[rec_idx] = mu * float(np.sum((lam ** (2 * alpha * np.arange(J))) * (u * u)))
            rec_idx += 1
        if step < nsteps:
            u = _rk4_step(u, dt, rhs_fn)
            if not np.all(np.isfinite(u)):
                return {
                    "t": t_hist[:rec_idx],
                    "u_hist": u_hist[:rec_idx],
                    "E_hist": E_hist[:rec_idx],
                    "sgs_power_hist": P_hist[:rec_idx],
                    "success": False,
                    "message": "NaN encountered",
                }
            t += dt

    return {
        "t": t_hist,
        "u_hist": u_hist,
        "E_hist": E_hist,
        "sgs_power_hist": P_hist,
        "success": True,
        "message": "OK",
    }
