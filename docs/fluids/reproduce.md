# Fluids evidence reproduction (commands)

## Fast checks (default CI-safe)
```
./check_python.sh
./check_lean.sh
python python/scripts/ns_run_all.py
```

## Regenerate full fluids evidence (medium)
```
NS_INCLUDE_SLOW=1 python python/scripts/ns_run_all.py
```

## Regenerate optional Track B evidence
```
NS_INCLUDE_OPTIONAL=1 python python/scripts/ns_run_all.py
```

## Full evidence (slow + optional)
```
NS_INCLUDE_SLOW=1 NS_INCLUDE_OPTIONAL=1 python python/scripts/ns_run_all.py
```

## Run slow tests explicitly
```
pytest -m slow -q
```

## Expected artifacts (fluids)
- ns05_energy_ledger.png
- ns06_uv_completion.png
- ns07_mu_convergence.png
- ns08_closure_kernel.png
- ns09_invariance.png
- ns14_energy_ledger_3d.png
- ns15_uv_completion_3d.png
- ns16_mu_convergence_3d.png
- ns18_exponent_sanity.png
- (slow) ns20_closure_kernel_3d.png
- (slow) ns21_closure_kernel_3d_forced_avg.png
- (optional) ns22_markov_laplacian_scaling.png
