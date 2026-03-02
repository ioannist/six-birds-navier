This file maps Six Primitives (P1–P6) to the BH interface-EFT pipeline.
It is a drafting aid; not final text.
Keep this compact and schematic.

| Primitive | One-line BH instantiation | Concrete repo anchor (file/script/theorem) | What the paper should claim |
|---|---|---|---|
| P1 Operator rewrite | Integrating out interior \rightarrow boundary operator \(\psi'(x_0,\omega)=Y(\omega)\psi(x_0,\omega)\). | docs/blackholes/conventions.md; python/src/bhwave/scattering.py | - Exterior dynamics close once Y(\(\omega\)) is specified. |
| P2 Constraints / admissible class | Constraints on Z(\(\omega\)) from passivity, causality, stability. | python/scripts/bh04_passivity_checks.py; python/scripts/bh05_kramers_kronig_numeric.py; python/scripts/bh10_cross_channel_linkage.py | - Most interior stories are ruled out as inadmissible boundary responses. |
| P3 Protocol holonomy | Inference/measurement protocol determines what is observable; two-point extraction vs naive ratio illustrates this. | python/scripts/bh16_fdtd_inference_demo.py; python/tests/test_bh16_two_point_estimator.py | - Protocol matters; correct observable extraction is part of the theory. |
| P4 Staging / sector decomposition | Regimes (absorbing vs reflective; non-superradiant vs superradiant active band) define sectors with different constraints. | python/scripts/bh11_superradiance_ledger_demo.py; python/scripts/bh14_boundary_stability_compare.py; python/scripts/bh17_rw_barrier_compare.py | - Different physical sectors correspond to different allowed regions in Z-space. |
| P5 Packaging / interface object | Interior complexity \rightarrow single interface response Z(\(\omega\)) (or R(\(\omega\))). | docs/blackholes/theorem_inventory.md; docs/blackholes/conventions.md; python/src/bhwave/scattering.py | - The interior is an equivalence class under the packaging map. |
| P6 Accounting / ledger | Flux ledger \(j=\omega(|R|^2-1)\) and causality/KK link; energy ledger in time-domain. | python/scripts/bh04_passivity_checks.py; formal/BlackHoles/Impedance.lean; python/scripts/bh13_energy_ledger_demo.py | - Energy/causality ledgers impose non-negotiable inequalities. |

Figure hooks (placeholders):
- python/artifacts/bh18_constraint_stacking.png
- python/artifacts/bh10_cross_channel.png
- python/artifacts/bh11_superradiance_gain.png
- python/artifacts/bh17_rw_vs_pt.png
