import numpy as np

from nswave.galerkin_core import (
    B_of_u,
    GalerkinModel,
    dissipation,
    energy,
    make_energy_preserving_tensor,
    make_psd_operator,
    simulate,
)


def _make_model(n: int, nu: float, mu: float, seed: int = 0) -> GalerkinModel:
    rng = np.random.default_rng(seed)
    T = rng.standard_normal((n, n, n))
    C = make_energy_preserving_tensor(T)
    A = make_psd_operator(n, seed=seed + 1)
    Aalpha = make_psd_operator(n, seed=seed + 2)
    return GalerkinModel(n=n, A=A, Aalpha=Aalpha, nu=nu, mu=mu, C=C)


def test_energy_cancellation_property():
    n = 12
    rng = np.random.default_rng(0)
    T = rng.standard_normal((n, n, n))
    C = make_energy_preserving_tensor(T)
    model = GalerkinModel(n=n, A=np.eye(n), Aalpha=np.eye(n), nu=0.0, mu=0.0, C=C)

    for _ in range(10):
        u = rng.standard_normal(n)
        Bu = B_of_u(model, u)
        val = float(u @ Bu)
        bound = 1e-12 * (np.linalg.norm(u) ** 3) + 1e-12
        assert abs(val) <= bound


def test_dissipation_nonnegative():
    n = 12
    model = _make_model(n, nu=0.5, mu=0.3, seed=0)
    rng = np.random.default_rng(1)
    for _ in range(10):
        u = rng.standard_normal(n)
        d = dissipation(model, u)
        assert d >= -1e-12


def test_energy_monotonicity_dissipative():
    n = 12
    model = _make_model(n, nu=0.5, mu=0.3, seed=0)
    rng = np.random.default_rng(2)
    u0 = rng.standard_normal(n)

    dt = 1e-3
    nsteps = 1000
    out = simulate(model, u0, dt, nsteps)
    E = out["E"]

    max_increase = np.max(E[1:] - E[:-1])
    assert max_increase <= 1e-6 * E[0] + 1e-10


def test_energy_conservative():
    n = 12
    model = _make_model(n, nu=0.0, mu=0.0, seed=0)
    rng = np.random.default_rng(3)
    u0 = rng.standard_normal(n)

    dt = 1e-3
    nsteps = 1000
    out = simulate(model, u0, dt, nsteps)
    E = out["E"]

    rel_err = np.max(np.abs(E - E[0])) / max(E[0], 1e-12)
    assert rel_err <= 1e-3
