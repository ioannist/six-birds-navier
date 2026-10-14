Draft figure manifest for BH section; captions are stubs.

- Figure ID: BH-F1
  Artifact filename: python/artifacts/bh05_kk_case2.png
  Supports claim: KK reconstruction of Re Z from Im Z for causal passive family.
  Intended location: sec:bh-crosschannel
  Caption stub: Re Z reconstructed from Im Z via KK on a Debye family. Shows agreement in the mid-band and clarifies truncation effects.
  Script: python/scripts/bh05_kramers_kronig_numeric.py

- Figure ID: BH-F2
  Artifact filename: python/artifacts/bh08_waveform.png
  Supports claim: Time-domain echoes from a realizable impedance boundary.
  Intended location: sec:bh-echoes
  Caption stub: Recorded waveform at an exterior observer for a cavity with a realizable boundary. Echo spacing is consistent with 2*Dx.
  Script: python/scripts/bh08_time_domain_demo.py

- Figure ID: BH-F3
  Artifact filename: python/artifacts/bh10_cross_channel.png
  Supports claim: Conservative/dissipative linkage via KK (cross-channel closure).
  Intended location: sec:bh-crosschannel
  Caption stub: Re Z predicted from Im Z on a causal Debye sum. Demonstrates the dispersion link between static response and dissipation.
  Script: python/scripts/bh10_cross_channel_linkage.py

- Figure ID: BH-F4
  Artifact filename: python/artifacts/bh11_superradiance_gain.png
  Supports claim: Toy superradiance ledger and gain > 1 band.
  Intended location: sec:bh-protocol
  Caption stub: Gain curve for a toy active-band reflectivity. Illustrates the ledger-aware instability hook when |R0 r|>1.
  Script: python/scripts/bh11_superradiance_ledger_demo.py

- Figure ID: BH-F5
  Artifact filename: python/artifacts/bh17_rw_vs_pt.png
  Supports claim: RW barrier realism vs P"oschl--Teller fit.
  Intended location: sec:bh-crosschannel
  Caption stub: |R| and |T| for RW (ell=2) compared to matched P"oschl--Teller. Confirms the toy barrier captures the peak neighborhood.
  Script: python/scripts/bh17_rw_barrier_compare.py

- Figure ID: BH-F6
  Artifact filename: python/artifacts/bh15_freqdep_echo_error.png
  Supports claim: Echo transfer formula holds for frequency-dependent Z(omega).
  Intended location: sec:bh-echoes
  Caption stub: Error between R_pred and ODE solver for Debye Z(omega). Validates the transfer formula in the dispersive case.
  Script: python/scripts/bh15_echo_freqdep_validation.py

- Figure ID: BH-F7
  Artifact filename: python/artifacts/bh13_energy.png
  Supports claim: Discrete energy ledger does not inject energy for passive boundary.
  Intended location: sec:bh-echoes
  Caption stub: Discrete energy time series for a passive Debye boundary. Energy is nonincreasing up to numerical tolerance.
  Script: python/scripts/bh13_energy_ledger_demo.py

- Figure ID: BH-F8
  Artifact filename: python/artifacts/bh16_fdtd_fit_waveform.png
  Supports claim: End-to-end inference recovers passive Z from FDTD waveform.
  Intended location: sec:bh-protocol
  Caption stub: True vs recovered waveform for the fitted Debye parameters. Shows time-domain agreement in the echo window.
  Script: python/scripts/bh16_fdtd_inference_demo.py

- Figure ID: BH-F9
  Artifact filename: python/artifacts/bh16_fdtd_fit_R_error.png
  Supports claim: Frequency-domain residuals for BH-16 inference.
  Intended location: sec:bh-protocol
  Caption stub: Spectral residual |R_meas - R_model| on the fit band. Confirms frequency-domain consistency of the recovered model.
  Script: python/scripts/bh16_fdtd_inference_demo.py

- Figure ID: BH-F10
  Artifact filename: python/artifacts/bh18_constraint_stacking.png
  Supports claim: Constraint stacking yields a feasible region in (a1,b1) space.
  Intended location: sec:bh-crosschannel
  Caption stub: Number of toy constraints satisfied over the Debye grid. Highlights the feasible region used for cross-channel closure.
  Script: python/scripts/bh18_constraint_stacking.py
