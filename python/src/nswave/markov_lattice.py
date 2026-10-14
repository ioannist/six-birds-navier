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
    if N <= 0 or M <= 0:
        raise ValueError("N and M must be positive")
    if not np.isfinite(rho0) or not 0 <= rho0 <= 1:
        raise ValueError("rho0 must lie in [0,1]")
    if not np.isfinite(eps) or abs(eps) > rho0:
        raise ValueError("|eps| must not exceed rho0")
    rng = np.random.default_rng(seed)
    x = np.arange(N)
    m = eps * np.cos(2 * np.pi * k_mode * x / N)
    p_plus = 0.5 * rho0 + 0.5 * m
    p_minus = 0.5 * rho0 - 0.5 * m
    p_empty = np.full_like(p_plus, 1.0 - rho0)

    probs = np.stack([p_empty, p_plus, p_minus], axis=0)

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
    """Apply two disjoint-bond matchings on an even periodic lattice.

    Each occupied/empty bond swaps at the species/direction probability.
    Fixed even-then-odd composition is the specified Markov update; it is
    not an exact continuous-time ASEP generator.
    """
    if state.ndim != 2 or state.shape[0] == 0 or state.shape[1] < 2 or state.shape[1] % 2:
        raise ValueError("state must be a nonempty ensemble on an even lattice N>=2")
    if np.any((state != 0) & (state != 1) & (state != -1)):
        raise ValueError("states must lie in {-1,0,1}")
    M, N = state.shape
    n_right = 0
    n_left = 0
    for probability in (p_plus, q_plus, p_minus, q_minus):
        arr = np.asarray(probability, dtype=float)
        if arr.shape not in ((), (N,), (M, N)):
            raise ValueError("probability shape must be scalar, (N,), or (M,N)")
        if not np.all(np.isfinite(arr)) or np.any(arr < 0) or np.any(arr > 1):
            raise ValueError("probabilities must be finite and lie in [0,1]")

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
    if snapshot_times is not None and list(snapshot_times):
        raise NotImplementedError("simulate returns only final state; record snapshots in the stepping loop")

    counters = {"n_right": 0, "n_left": 0}
    for step in range(n_steps):
        n_r, n_l = step_parallel_two_species_asep(
            state, p_plus, q_plus, p_minus, q_minus, rng
        )
        counters["n_right"] += n_r
        counters["n_left"] += n_l

    return state, counters
