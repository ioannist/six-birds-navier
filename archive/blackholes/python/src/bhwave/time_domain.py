"""Time-domain FDTD solver with a realizable impedance boundary."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np


@dataclass
class ImpedanceBoundary:
    a0: float
    a: np.ndarray
    b: np.ndarray
    q: np.ndarray

    @classmethod
    def from_constant_Z(cls, Z: float) -> "ImpedanceBoundary":
        return cls(a0=float(Z), a=np.zeros(0), b=np.zeros(0), q=np.zeros(0))

    @classmethod
    def from_debye(cls, a0: float, a: list[float], b: list[float]) -> "ImpedanceBoundary":
        a_arr = np.asarray(a, dtype=float)
        b_arr = np.asarray(b, dtype=float)
        if a_arr.shape != b_arr.shape:
            raise ValueError("a and b must have the same length")
        if np.any(a_arr < 0):
            raise ValueError("a must be nonnegative")
        if np.any(b_arr <= 0):
            raise ValueError("b must be positive")
        q = np.zeros_like(a_arr)
        return cls(a0=float(a0), a=a_arr, b=b_arr, q=q)

    def derivative_target(self, psi_t: float) -> float:
        return float(self.a0 * psi_t + np.sum(self.a * self.q))

    def step(self, psi_t: float, dt: float) -> None:
        if self.q.size == 0:
            return
        self.q = self.q + dt * (-self.b * self.q + self.b * psi_t)

    def step_exact(self, psi_t: float, dt: float) -> None:
        if self.q.size == 0:
            return
        decay = np.exp(-self.b * dt)
        self.q = decay * self.q + (1.0 - decay) * psi_t


@dataclass
class TimeDomainResult:
    x: np.ndarray
    t: np.ndarray
    psi_obs: np.ndarray
    x_obs: float
    dx: float
    dt: float
    x0: float
    x_max: float
    notes: str
    psi_full: np.ndarray | None = None
    t_full: np.ndarray | None = None
    q_full: np.ndarray | None = None
    x_probe_list: np.ndarray | None = None
    psi_probe: np.ndarray | None = None


def _laplacian_with_ghost(psi: np.ndarray, ghost_left: float, ghost_right: float, dx: float) -> np.ndarray:
    lap = np.empty_like(psi)
    inv_dx2 = 1.0 / (dx * dx)
    lap[1:-1] = (psi[2:] - 2.0 * psi[1:-1] + psi[:-2]) * inv_dx2
    lap[0] = (psi[1] - 2.0 * psi[0] + ghost_left) * inv_dx2
    lap[-1] = (ghost_right - 2.0 * psi[-1] + psi[-2]) * inv_dx2
    return lap


def discrete_energy(
    psi_curr: np.ndarray,
    psi_prev: np.ndarray,
    V: np.ndarray,
    dx: float,
    dt: float,
    *,
    boundary: ImpedanceBoundary | None = None,
    gamma: np.ndarray | None = None,
) -> float:
    psi_curr = np.asarray(psi_curr, dtype=float)
    psi_prev = np.asarray(psi_prev, dtype=float)
    V = np.asarray(V, dtype=float)
    if psi_curr.shape != psi_prev.shape or psi_curr.shape != V.shape:
        raise ValueError("psi_curr, psi_prev, and V must have the same shape")
    if dx <= 0 or dt <= 0:
        raise ValueError("dx and dt must be positive")

    psi_t = (psi_curr - psi_prev) / dt
    psi_x = np.empty_like(psi_curr)
    psi_x[1:-1] = (psi_curr[2:] - psi_curr[:-2]) / (2.0 * dx)
    psi_x[0] = (psi_curr[1] - psi_curr[0]) / dx
    psi_x[-1] = (psi_curr[-1] - psi_curr[-2]) / dx

    energy_density = 0.5 * (psi_t**2 + psi_x**2 + V * psi_curr**2)
    energy = float(np.sum(energy_density) * dx)

    if boundary is not None and boundary.q.size > 0:
        energy += float(0.5 * np.sum((boundary.a / boundary.b) * (boundary.q**2)))

    return energy


def simulate_wave_fdtd(
    *,
    V_func: Callable[[np.ndarray], np.ndarray],
    x0: float,
    x_max: float,
    dx: float,
    dt: float,
    t_max: float,
    boundary: ImpedanceBoundary,
    x_obs: float,
    sponge_start: float | None = None,
    sponge_gamma_max: float = 0.0,
    sponge_power: float = 2.0,
    psi_init: Callable[[np.ndarray], np.ndarray] | None = None,
    psi_t_init: Callable[[np.ndarray], np.ndarray] | None = None,
    record_every: int = 1,
    record_full_every: int | None = None,
    x_probe_list: list[float] | None = None,
    boundary_update: str = "exact",
    boundary_bc: str = "simple",
) -> TimeDomainResult:
    if dt <= 0 or dx <= 0:
        raise ValueError("dt and dx must be positive")
    if dt / dx > 1.0:
        raise ValueError("CFL violated: dt/dx must be <= 1.0")
    if x_max <= x0:
        raise ValueError("x_max must be greater than x0")
    if record_every <= 0:
        raise ValueError("record_every must be positive")
    if record_full_every is not None and record_full_every <= 0:
        raise ValueError("record_full_every must be positive")
    if boundary_update not in {"exact", "euler"}:
        raise ValueError("boundary_update must be 'exact' or 'euler'")
    if boundary_bc not in {"simple", "predictor_corrector"}:
        raise ValueError("boundary_bc must be 'simple' or 'predictor_corrector'")

    x = np.arange(x0, x_max + 0.5 * dx, dx, dtype=float)
    if not (x0 <= x_obs <= x_max):
        raise ValueError("x_obs must be within [x0, x_max]")
    j_obs = int(np.argmin(np.abs(x - x_obs)))

    probe_list = None
    probe_indices = None
    if x_probe_list is not None:
        probe_list = np.asarray(x_probe_list, dtype=float)
        if np.any((probe_list < x0) | (probe_list > x_max)):
            raise ValueError("x_probe_list entries must be within [x0, x_max]")
        probe_indices = np.array([int(np.argmin(np.abs(x - xp))) for xp in probe_list], dtype=int)

    V = np.asarray(V_func(x), dtype=float)
    if V.shape != x.shape:
        raise ValueError("V_func must return an array of the same shape as x")

    gamma = np.zeros_like(x)
    if sponge_start is not None and sponge_gamma_max > 0.0:
        if sponge_start >= x_max:
            raise ValueError("sponge_start must be less than x_max")
        mask = x > sponge_start
        if np.any(mask):
            scale = (x[mask] - sponge_start) / (x_max - sponge_start)
            gamma[mask] = sponge_gamma_max * (scale**sponge_power)

    if psi_init is None:
        psi0 = np.zeros_like(x)
    else:
        psi0 = np.asarray(psi_init(x), dtype=float)
    if psi_t_init is None:
        v0 = np.zeros_like(x)
    else:
        v0 = np.asarray(psi_t_init(x), dtype=float)
    if psi0.shape != x.shape or v0.shape != x.shape:
        raise ValueError("psi_init and psi_t_init must match the x grid shape")

    psi_t0 = float(v0[0])
    dtarget = boundary.derivative_target(psi_t0)
    ghost_left = psi0[1] - 2.0 * dx * dtarget
    psi_tR = float(v0[-1])
    ghost_right = psi0[-2] - 2.0 * dx * psi_tR
    lap0 = _laplacian_with_ghost(psi0, ghost_left, ghost_right, dx)
    psi1 = psi0 + dt * v0 + 0.5 * dt * dt * (lap0 - V * psi0 - gamma * v0)

    psi_prev = psi0
    psi_curr = psi1

    n_steps = int(np.round(t_max / dt))
    t_rec: list[float] = []
    psi_rec: list[float] = []
    psi_probe: list[list[float]] | None = None
    if probe_indices is not None:
        psi_probe = [[] for _ in range(len(probe_indices))]
    psi_full: list[np.ndarray] | None = [] if record_full_every is not None else None
    t_full: list[float] | None = [] if record_full_every is not None else None
    q_full: list[np.ndarray] | None = [] if record_full_every is not None else None

    def record(step_idx: int, psi: np.ndarray) -> None:
        if step_idx % record_every == 0:
            t_rec.append(step_idx * dt)
            psi_rec.append(float(psi[j_obs]))
            if psi_probe is not None and probe_indices is not None:
                for k, idx in enumerate(probe_indices):
                    psi_probe[k].append(float(psi[idx]))

    def record_full(step_idx: int, psi: np.ndarray) -> None:
        if psi_full is None or t_full is None or q_full is None:
            return
        if step_idx % record_full_every == 0:
            t_full.append(step_idx * dt)
            psi_full.append(psi.copy())
            q_full.append(boundary.q.copy())

    record(0, psi_prev)
    record(1, psi_curr)
    record_full(0, psi_prev)
    record_full(1, psi_curr)

    for step in range(1, n_steps):
        psi_t0_pred = float((psi_curr[0] - psi_prev[0]) / dt)
        if boundary_update == "exact" and boundary.q.size > 0:
            decay = np.exp(-boundary.b * dt)
            q_next_pred = decay * boundary.q + (1.0 - decay) * psi_t0_pred
            dtarget_pred = float(boundary.a0 * psi_t0_pred + np.sum(boundary.a * q_next_pred))
        else:
            dtarget_pred = boundary.derivative_target(psi_t0_pred)
        ghost_left = psi_curr[1] - 2.0 * dx * dtarget_pred
        psi_tR = float((psi_curr[-1] - psi_prev[-1]) / dt)
        ghost_right = psi_curr[-2] - 2.0 * dx * psi_tR

        lap = _laplacian_with_ghost(psi_curr, ghost_left, ghost_right, dx)
        psi_next = (
            2.0 * psi_curr
            - psi_prev
            + (dt * dt) * lap
            - (dt * dt) * V * psi_curr
            - gamma * dt * (psi_curr - psi_prev)
        )
        sommerfeld_coef = (dt - dx) / (dt + dx)
        psi_next[-1] = psi_curr[-2] + sommerfeld_coef * (psi_next[-2] - psi_curr[-1])

        psi_t0_used = psi_t0_pred
        if boundary_bc == "predictor_corrector":
            if boundary_update == "exact":
                A = boundary.a0
                B = 0.0
                if boundary.q.size > 0:
                    decay = np.exp(-boundary.b * dt)
                    A += float(np.sum(boundary.a * (1.0 - decay)))
                    B = float(np.sum(boundary.a * (decay * boundary.q)))
                lap_base = 2.0 * (psi_curr[1] - psi_curr[0]) / (dx * dx)
                const = (
                    2.0 * psi_curr[0]
                    - psi_prev[0]
                    + (dt * dt) * (lap_base - (2.0 / dx) * B)
                    - (dt * dt) * V[0] * psi_curr[0]
                    - gamma[0] * dt * (psi_curr[0] - psi_prev[0])
                )
                alpha = (dt / dx) * A
                psi_next0 = (const + alpha * psi_prev[0]) / (1.0 + alpha)
                psi_next[0] = psi_next0
                psi_t0_corr = float((psi_next0 - psi_prev[0]) / (2.0 * dt))
                if boundary.q.size > 0:
                    q_next_corr = decay * boundary.q + (1.0 - decay) * psi_t0_corr
                psi_t0_used = psi_t0_corr
            else:
                psi_t0_corr = float((psi_next[0] - psi_prev[0]) / (2.0 * dt))
                dtarget_corr = boundary.derivative_target(psi_t0_corr)
                ghost_left_corr = psi_curr[1] - 2.0 * dx * dtarget_corr
                lap0 = (psi_curr[1] - 2.0 * psi_curr[0] + ghost_left_corr) / (dx * dx)
                psi_next[0] = (
                    2.0 * psi_curr[0]
                    - psi_prev[0]
                    + (dt * dt) * lap0
                    - (dt * dt) * V[0] * psi_curr[0]
                    - gamma[0] * dt * (psi_curr[0] - psi_prev[0])
                )
                psi_t0_used = psi_t0_corr

        if boundary_update == "exact":
            if boundary.q.size > 0:
                if boundary_bc == "predictor_corrector" and "q_next_corr" in locals():
                    boundary.q = q_next_corr
                elif "q_next_pred" in locals():
                    boundary.q = q_next_pred
                else:
                    boundary.step_exact(psi_t0_used, dt)
        else:
            boundary.step(psi_t0_used, dt)

        psi_prev, psi_curr = psi_curr, psi_next
        record(step + 1, psi_curr)
        record_full(step + 1, psi_curr)

    return TimeDomainResult(
        x=x,
        t=np.asarray(t_rec, dtype=float),
        psi_obs=np.asarray(psi_rec, dtype=float),
        x_obs=float(x_obs),
        dx=float(dx),
        dt=float(dt),
        x0=float(x0),
        x_max=float(x_max),
        notes="FDTD with impedance and Sommerfeld boundaries",
        psi_full=None if psi_full is None else np.asarray(psi_full, dtype=float),
        t_full=None if t_full is None else np.asarray(t_full, dtype=float),
        q_full=None if q_full is None else np.asarray(q_full, dtype=float),
        x_probe_list=None if probe_list is None else probe_list.copy(),
        psi_probe=None if psi_probe is None else np.asarray(psi_probe, dtype=float),
    )
