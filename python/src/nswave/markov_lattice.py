"""Simple lattice Markov substrate for Track B (two-species ASEP with exclusion)."""

from __future__ import annotations

from typing import Iterable

import numpy as np

Array = np.ndarray


def init_sinusoidal_imbalance_ensemble(
    N: int,
    M: int,
    k_mode: int,
    rho0: float,
    eps: float,
    seed: int = 0,
) -> Array:
    """Initialize an ensemble with sinusoidal momentum imbalance."""
    rng = np.random.default_rng(seed)
    x = np.arange(N)
    m = eps * np.cos(2 * np.pi * k_mode * x / N)
    p_plus = 0.5 * rho0 + 0.5 * m
    p_minus = 0.5 * rho0 - 0.5 * m
    p_empty = np.full_like(p_plus, 1.0 - rho0)

    p_plus = np.clip(p_plus, 0.0, 1.0)
    p_minus = np.clip(p_minus, 0.0, 1.0)
    p_empty = np.clip(p_empty, 0.0, 1.0)

    probs = np.stack([p_empty, p_plus, p_minus], axis=0)
    probs /= np.sum(probs, axis=0, keepdims=True)

    state = np.zeros((M, N), dtype=np.int8)
    draws = rng.random((M, N))
    cdf0 = probs[0]
    cdf1 = probs[0] + probs[1]
    state[draws > cdf0] = 1
    state[draws > cdf1] = -1
    return state


def _apply_moves(
    state: Array,
    i_idx: Array,
    j_idx: Array,
    move_mask: Array,
    value: int,
) -> int:
    if not np.any(move_mask):
        return 0
    m_idx, b_idx = np.where(move_mask)
    state[m_idx, i_idx[b_idx]] = 0
    state[m_idx, j_idx[b_idx]] = value
    return int(move_mask.sum())


def step_parallel_two_species_asep(
    state: Array,
    p_plus: float,
    q_plus: float,
    p_minus: float,
    q_minus: float,
    rng: np.random.Generator,
) -> tuple[int, int]:
    """Advance the lattice by one step using even/odd bond updates."""
    M, N = state.shape
    n_right = 0
    n_left = 0

    def _prob_block(prob: float | Array, idx: Array) -> Array | float:
        prob_arr = np.asarray(prob)
        if prob_arr.ndim == 0:
            return prob
        if prob_arr.ndim == 1:
            if prob_arr.shape[0] != N:
                raise ValueError("1D probability arrays must have length N")
            return np.broadcast_to(prob_arr[idx], (M, idx.size))
        if prob_arr.ndim == 2:
            return prob_arr[:, idx]
        raise ValueError("probability array must be scalar, (N,), or (M,N)")

    for start in (0, 1):
        indices = np.arange(start, N, 2)
        j_idx = (indices + 1) % N
        s_i = state[:, indices]
        s_j = state[:, j_idx]

        empty_i = s_i == 0
        empty_j = s_j == 0

        plus_i = s_i == 1
        minus_i = s_i == -1
        plus_j = s_j == 1
        minus_j = s_j == -1

        p_plus_block = _prob_block(p_plus, indices)
        q_plus_block = _prob_block(q_plus, j_idx)
        p_minus_block = _prob_block(p_minus, indices)
        q_minus_block = _prob_block(q_minus, j_idx)

        # Right moves from i -> j.
        r_plus = plus_i & empty_j & (rng.random(s_i.shape) < p_plus_block)
        r_minus = minus_i & empty_j & (rng.random(s_i.shape) < p_minus_block)

        # Left moves from j -> i.
        l_plus = plus_j & empty_i & (rng.random(s_j.shape) < q_plus_block)
        l_minus = minus_j & empty_i & (rng.random(s_j.shape) < q_minus_block)

        n_right += _apply_moves(state, indices, j_idx, r_plus, 1)
        n_right += _apply_moves(state, indices, j_idx, r_minus, -1)

        n_left += _apply_moves(state, j_idx, indices, l_plus, 1)
        n_left += _apply_moves(state, j_idx, indices, l_minus, -1)

    return n_right, n_left


def simulate(
    state: Array,
    *,
    p_plus: float,
    q_plus: float,
    p_minus: float,
    q_minus: float,
    n_steps: int,
    rng: np.random.Generator,
    snapshot_times: Iterable[int] | None = None,
) -> tuple[Array, dict[str, int]]:
    """Run the Markov chain for n_steps and return the final state."""
    if n_steps < 0:
        raise ValueError("n_steps must be nonnegative")
    snapshot_set = set(snapshot_times or [])

    counters = {"n_right": 0, "n_left": 0}
    for step in range(n_steps):
        n_r, n_l = step_parallel_two_species_asep(
            state, p_plus, q_plus, p_minus, q_minus, rng
        )
        counters["n_right"] += n_r
        counters["n_left"] += n_l
        if step in snapshot_set:
            pass

    return state, counters
