"""Finite-dimensional Galerkin core with NS-style energy bookkeeping."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np


@dataclass
class GalerkinModel:
    n: int
    A: np.ndarray
    Aalpha: np.ndarray
    nu: float
    mu: float
    C: np.ndarray


def make_energy_preserving_tensor(T: np.ndarray) -> np.ndarray:
    """Skew-symmetrize the first two indices to enforce u·B(u,u)=0."""
    return T - np.swapaxes(T, 0, 1)


def B_of_u(model: GalerkinModel, u: np.ndarray) -> np.ndarray:
    """Quadratic nonlinearity B(u,u) defined by the tensor C."""
    return np.einsum("ijk,j,k->i", model.C, u, u)


def rhs(model: GalerkinModel, u: np.ndarray, f: Optional[np.ndarray] = None) -> np.ndarray:
    """Right-hand side for the Galerkin ODE."""
    nonlinear = B_of_u(model, u)
    linear = model.nu * (model.A @ u) + model.mu * (model.Aalpha @ u)
    if f is None:
        return -nonlinear - linear
    return -nonlinear - linear + f


def energy(u: np.ndarray) -> float:
    """Kinetic energy 0.5*||u||^2."""
    return 0.5 * float(u @ u)


def dissipation(model: GalerkinModel, u: np.ndarray) -> float:
    """Dissipation rate induced by A and Aalpha."""
    return model.nu * float(u @ (model.A @ u)) + model.mu * float(u @ (model.Aalpha @ u))


def step_rk4(
    model: GalerkinModel, u: np.ndarray, dt: float, f: Optional[np.ndarray] = None
) -> np.ndarray:
    """One RK4 step for the Galerkin ODE."""
    k1 = rhs(model, u, f)
    k2 = rhs(model, u + 0.5 * dt * k1, f)
    k3 = rhs(model, u + 0.5 * dt * k2, f)
    k4 = rhs(model, u + dt * k3, f)
    return u + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def simulate(
    model: GalerkinModel,
    u0: np.ndarray,
    dt: float,
    nsteps: int,
    f: Optional[np.ndarray] = None,
) -> dict:
    """Simulate and return time, energy, and dissipation histories."""
    u = u0.copy()
    t = np.zeros(nsteps + 1)
    E = np.zeros(nsteps + 1)
    D = np.zeros(nsteps + 1)
    E[0] = energy(u)
    D[0] = dissipation(model, u)

    for i in range(nsteps):
        u = step_rk4(model, u, dt, f)
        t[i + 1] = (i + 1) * dt
        E[i + 1] = energy(u)
        D[i + 1] = dissipation(model, u)

    return {"t": t, "E": E, "D": D}


def make_psd_operator(
    n: int, seed: int = 0, *, rank: Optional[int] = None, scale: float = 1.0
) -> np.ndarray:
    """Construct a symmetric PSD operator via M M^T."""
    rng = np.random.default_rng(seed)
    r = n if rank is None else int(rank)
    if r <= 0:
        raise ValueError("rank must be positive")

    M = rng.standard_normal((n, r))
    A = (M @ M.T) / float(r)
    A = 0.5 * (A + A.T)
    return scale * A
