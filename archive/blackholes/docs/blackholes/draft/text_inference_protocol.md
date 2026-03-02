# Inference as protocol (P3): extracting Z(ω) from signals

## Why protocol matters
- TODO(BH): “observable” depends on measurement/extraction map, not just dynamics
- TODO(BH): explain free-region plane-wave decomposition as a protocol choice

## BH-16 pipeline: waveform → R_out(ω) → fit Debye Z(ω) → replay
- TODO(BH): two-point extraction; define y=x-x0; solve for A_in,A_out; R_out=A_out/A_in
- TODO(BH): passive family fit: Z(ω)=a0+a1/(1-iω/b1), enforce a0,a1,b1>0
- TODO(BH): replay FDTD with recovered params; compare waveform in echo window

## Failure mode: why the baseline ratio estimator can be ill-posed
- TODO(BH): explain factorization assumption F≈S·R; when it fails (probe doesn’t see incident; source differs)
- TODO(BH): mention we fixed by moving x_src and using two-point estimator + band restriction

## What we claim vs what we do not
- TODO(BH): claim: feasibility of recovering passive Z in toy setting; not claim: observational inference for real BH

## Evidence pointers
- TODO(BH): cite bh16_fdtd_inference_demo.py, inference.py helper, tests, artifacts
