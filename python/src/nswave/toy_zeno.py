from __future__ import annotations

import numpy as np


def toy_case_series(case_id: str, J: int) -> dict[str, np.ndarray]:
    """
    Return toy series for Zeno/No-Zeno cases.

    Keys: j, w, Cap, delta_t, t_J, sum_w_over_cap_J.
    """
    if J <= 0:
        raise ValueError("J must be positive")

    j = np.arange(J, dtype=float)
    case_id_norm = case_id.strip()

    if case_id_norm == "Z1_fast_capacity":
        w = np.ones(J, dtype=float)
        cap = 2.0 ** j
    elif case_id_norm == "Z2_vanishing_work":
        w = 2.0 ** (-j)
        cap = np.ones(J, dtype=float)
    elif case_id_norm == "NZ_polynomial_capacity":
        w = np.ones(J, dtype=float)
        cap = j + 1.0
    else:
        raise ValueError(f"Unknown case_id: {case_id}")

    delta_t = w / cap
    t_J = np.cumsum(delta_t)
    sum_w_over_cap_J = np.cumsum(w / cap)

    return {
        "j": j,
        "w": w,
        "Cap": cap,
        "delta_t": delta_t,
        "t_J": t_J,
        "sum_w_over_cap_J": sum_w_over_cap_J,
    }
