#!/usr/bin/env python3
from __future__ import annotations

import pathlib
import sys

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

from bhwave.barrier import barrier_reflection_transmission, poeschl_teller_V  # noqa: E402
from bhwave.inference import (  # noqa: E402
    R_pred_grid,
    Z_debye_k1,
    estimate_Rout_from_two_point_spectra,
    fit_debye_k1_from_Rpred,
    passivity_residuals,
)
from bhwave.time_domain import ImpedanceBoundary, simulate_wave_fdtd  # noqa: E402


def run_fdtd(
    *,
    V_func,
    x0: float,
    x_max: float,
    dx: float,
    dt: float,
    t_max: float,
    x_obs: float,
    sponge_start: float | None,
    sponge_gamma_max: float,
    boundary: ImpedanceBoundary,
    boundary_update: str,
    boundary_bc: str,
    psi_init,
    psi_t_init,
    x_probe_list: list[float] | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    result = simulate_wave_fdtd(
        V_func=V_func,
        x0=x0,
        x_max=x_max,
        dx=dx,
        dt=dt,
        t_max=t_max,
        boundary=boundary,
        x_obs=x_obs,
        sponge_start=sponge_start,
        sponge_gamma_max=sponge_gamma_max,
        psi_init=psi_init,
        psi_t_init=psi_t_init,
        record_every=1,
        x_probe_list=x_probe_list,
        boundary_update=boundary_update,
        boundary_bc=boundary_bc,
    )
    return result.t, result.psi_obs, result.psi_probe


def main() -> None:
    V0 = 1.0
    a = 1.0
    x_ref = 0.0
    x_peak = 0.0

    def V_func(x: np.ndarray) -> np.ndarray:
        return V0 / np.cosh(a * (x - x_peak)) ** 2

    def V_scalar(x: float) -> float:
        return poeschl_teller_V(x, V0, a, x_peak)

    x0 = -8.0
    x_max = 40.0
    x_obs = 20.0
    Dx = x_ref - x0
    sponge_start = 30.0
    sponge_gamma_max = 1.0

    dx = 0.02
    dt = 0.012
    t_max = 70.0

    boundary_update = "exact"
    boundary_bc = "predictor_corrector"

    x_src = 28.0
    sigma = 0.8
    k0 = 1.5

    def psi_init(x: np.ndarray) -> np.ndarray:
        gaussian = np.exp(-((x - x_src) ** 2) / (2.0 * sigma**2))
        return gaussian * np.cos(k0 * (x - x_src))

    def psi_t_init(x: np.ndarray) -> np.ndarray:
        gaussian = np.exp(-((x - x_src) ** 2) / (2.0 * sigma**2))
        d_gauss = -((x - x_src) / (sigma**2)) * gaussian
        phase = k0 * (x - x_src)
        dpsi_dx = d_gauss * np.cos(phase) - k0 * gaussian * np.sin(phase)
        # Left-moving packet for u_tt = u_xx uses psi_t = +psi_x.
        return dpsi_dx

    a0_true, a1_true, b1_true = 0.2, 0.5, 2.0
    p0 = (0.15, 0.3, 1.0)

    x2 = x_obs + dx
    probe_list = [x_obs, x2]

    print("Config:")
    print(f"  dx={dx}, dt={dt}, t_max={t_max}, Dx={Dx}")
    print(f"  x_src={x_src}, x_obs={x_obs}, x2={x2}, k0={k0}")
    print(f"  boundary_update={boundary_update}, boundary_bc={boundary_bc}")
    print(f"  true params: a0={a0_true}, a1={a1_true}, b1={b1_true}")
    print(f"  init guess: {p0}")
    print()

    boundary_true = ImpedanceBoundary.from_debye(a0=a0_true, a=[a1_true], b=[b1_true])
    t_true, y_true, psi_probe_true = run_fdtd(
        V_func=V_func,
        x0=x0,
        x_max=x_max,
        dx=dx,
        dt=dt,
        t_max=t_max,
        x_obs=x_obs,
        sponge_start=sponge_start,
        sponge_gamma_max=sponge_gamma_max,
        boundary=boundary_true,
        boundary_update=boundary_update,
        boundary_bc=boundary_bc,
        psi_init=psi_init,
        psi_t_init=psi_t_init,
        x_probe_list=probe_list,
    )

    if psi_probe_true is None or psi_probe_true.shape[0] < 2:
        raise RuntimeError("probe recording failed")

    N = len(y_true)
    omega_pos = 2.0 * np.pi * np.fft.rfftfreq(N, dt)
    # rfft returns coefficients for an inverse series with exp(+i omega t);
    # conjugate to interpret amplitudes under the exp(-i omega t) convention.
    Psi1 = np.conjugate(np.fft.rfft(psi_probe_true[0]))
    Psi2 = np.conjugate(np.fft.rfft(psi_probe_true[1]))

    y1 = x_obs - x0
    y2 = x2 - x0
    R_meas_all, mask_amp = estimate_Rout_from_two_point_spectra(
        omega_pos=omega_pos,
        Psi1=Psi1,
        Psi2=Psi2,
        y1=y1,
        y2=y2,
        min_amp_frac=2e-2,
    )
    mask_band = (omega_pos >= 0.9) & (omega_pos <= 3.0)
    mask_passive = np.abs(R_meas_all) <= 1.2
    mask = mask_amp & mask_band & mask_passive
    if not np.any(mask):
        raise RuntimeError("no usable FFT bins for fitting")

    idx = np.where(mask)[0]
    max_points = 120
    if idx.size > max_points:
        sel = np.linspace(0, idx.size - 1, max_points).astype(int)
        idx_fit = idx[sel]
    else:
        idx_fit = idx

    omega_fit = omega_pos[idx_fit]
    R_meas = R_meas_all[idx_fit]

    r_fit = np.empty_like(omega_fit, dtype=complex)
    t_fit = np.empty_like(omega_fit, dtype=complex)
    for i, omega in enumerate(omega_fit):
        barrier = barrier_reflection_transmission(
            float(omega),
            V_scalar,
            x_min=-12.0,
            x_max=12.0,
            x_ref=x_ref,
        )
        r_fit[i] = barrier.R_barrier
        t_fit[i] = barrier.T_barrier

    fit = fit_debye_k1_from_Rpred(
        omega_pos=omega_fit,
        Rpred_target=R_meas,
        r=r_fit,
        t=t_fit,
        Dx=Dx,
        p0=p0,
        weight=np.ones_like(omega_fit),
    )

    a0_fit = fit["a0"]
    a1_fit = fit["a1"]
    b1_fit = fit["b1"]

    rel_err_a0 = abs(a0_fit - a0_true) / a0_true
    rel_err_a1 = abs(a1_fit - a1_true) / a1_true
    rel_err_b1 = abs(b1_fit - b1_true) / b1_true

    Z_fit = Z_debye_k1(omega_fit, a0_fit, a1_fit, b1_fit)
    R0_fit = (1.0 - Z_fit) / (1.0 + Z_fit)
    R_model_fit = R_pred_grid(omega_fit, R0_fit, r_fit, t_fit, Dx)

    resid_R, resid_Z = passivity_residuals(omega_fit, Z_fit, R0_fit)
    max_abs_err_R = float(np.max(np.abs(R_meas - R_model_fit)))

    print("Fit results:")
    print(f"  recovered a0={a0_fit:.6g} (rel err {rel_err_a0:.3%})")
    print(f"  recovered a1={a1_fit:.6g} (rel err {rel_err_a1:.3%})")
    print(f"  recovered b1={b1_fit:.6g} (rel err {rel_err_b1:.3%})")
    print(f"  passivity residuals: resid_R={resid_R:.3e}, resid_Z={resid_Z:.3e}")
    print(f"  max_abs_err_R={max_abs_err_R:.6g}")
    print()

    boundary_fit = ImpedanceBoundary.from_debye(a0=a0_fit, a=[a1_fit], b=[b1_fit])
    t_fit_wave, y_fit, _ = run_fdtd(
        V_func=V_func,
        x0=x0,
        x_max=x_max,
        dx=dx,
        dt=dt,
        t_max=t_max,
        x_obs=x_obs,
        sponge_start=sponge_start,
        sponge_gamma_max=sponge_gamma_max,
        boundary=boundary_fit,
        boundary_update=boundary_update,
        boundary_bc=boundary_bc,
        psi_init=psi_init,
        psi_t_init=psi_t_init,
        x_probe_list=probe_list,
    )

    if t_fit_wave.shape != t_true.shape:
        raise RuntimeError("fit waveform time grid mismatch")

    mask_window = t_true >= 20.0
    y_true_win = y_true[mask_window]
    y_fit_win = y_fit[mask_window]
    rel_L2 = np.linalg.norm(y_fit_win - y_true_win) / max(1e-12, np.linalg.norm(y_true_win))
    print(f"rel_L2(t>=20)={rel_L2:.6g}")

    artifacts_dir = ROOT / "artifacts"
    artifacts_dir.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(6.0, 4.0))
    plt.plot(t_true, y_true, label="true")
    plt.plot(t_fit_wave, y_fit, label="fit", alpha=0.8)
    plt.xlabel("t")
    plt.ylabel("psi_obs")
    plt.title("BH-16 FDTD waveform fit")
    plt.legend()
    plt.tight_layout()
    waveform_path = artifacts_dir / "bh16_fdtd_fit_waveform.png"
    plt.savefig(waveform_path)

    plt.figure(figsize=(6.0, 4.0))
    plt.plot(omega_fit, np.abs(R_meas - R_model_fit))
    plt.xlabel("omega")
    plt.ylabel("|R_meas - R_model|")
    plt.title("BH-16 frequency-domain fit residual")
    plt.tight_layout()
    error_path = artifacts_dir / "bh16_fdtd_fit_R_error.png"
    plt.savefig(error_path)

    print(f"Saved plot: {waveform_path}")
    print(f"Saved plot: {error_path}")

    if (
        rel_err_a0 > 0.2
        or rel_err_a1 > 0.2
        or rel_err_b1 > 0.2
        or max_abs_err_R > 5e-2
    ):
        print("FAIL: parameter recovery or spectral error out of tolerance")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
