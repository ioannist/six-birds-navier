# Six Birds: Navier--Stokes Instantiation

This repository contains the **Navier--Stokes instantiation** for the paper:

> **Six Birds for Navier--Stokes: A Mechanized No-Zeno Scaffold**
>
> Archived at: TBD
>
> DOI: TBD

This paper is the Navier--Stokes-focused instantiation of the emergence/interface calculus introduced in *Six Birds: Foundations of Emergence Calculus*. It formalizes a non-claiming decision scaffold: mechanized abstract inequalities in Lean, plus reproducible diagnostics and certificates in Python.

## Current conditional Navier theorem

The later conditional Navier manuscript is maintained in the sibling
`six-birds-needles` repository at `paper/needles_ns/main.tex`. Its current
spectral Hilbert XI and physical PDE Lean sources are under
`lean/SixBirdsNeedles/NSCore/`, principally `SharedXiStrongClosure.lean` and
`PhysicalClassicalReturn.lean`. That paper's older physical-window/BG route
has a separate hypothesis. This repository's No-Zeno scaffold and its
`review/implementation/needles/` copy are historical campaign artifacts,
not the source of the current conditional theorem.

## What this repository provides

The Navier--Stokes instantiation implements:

- **Lean/formal anchors**: machine-checked theorem-chain components (No-Zeno increment logic, work/capacity algebra, route-capacity polynomial bounds, finite-channel anti-localization lemma, and ICAP no-go under localization-rich tests)
- **Manifest-driven figure pipeline**: deterministic paper artifacts from `paper/figures_manifest.json` with paired `.png` + `.npz` outputs
- **SBT-legality diagnostics**: operational anti-localization/channel metrics and sweep scripts for stress-testing the feasibility hinge
- **Illegal-vs-legal certificate demo**: explicit localized-shell constructions that fail the same certificate passed by baseline initial data
- **Appendix autogen contract**: paper appendix figure table generated from manifest (`python/scripts/ns_paper_appendix_figure_table.py`)

## Scope and limitations

The paper is explicit about what it does and does not establish:

- It does **not** claim a proof of Clay 3D Navier--Stokes regularity
- Mechanized results are abstract/discrete inequality logic, not PDE closure theorems
- Numerical outputs are diagnostics and falsification hooks, not worst-case proofs
- The key NS-facing hinge (feasibility/channelization/anti-localization for canonical packaged bridge objects) remains open

## Install

```bash
python -m pip install -r python/requirements.txt
cd formal && lake build
```

## Test

```bash
./check_python.sh
./check_lean.sh
./check_paper.sh
```

## Run experiments

```bash
python python/scripts/ns_run_all.py
NS_INCLUDE_SLOW=1 python python/scripts/ns_run_all.py
NS_INCLUDE_OPTIONAL=1 python python/scripts/ns_run_all.py
python python/scripts/ns_paper_figures.py --include-slow --clean
```

## Build paper

```bash
bash paper/build_paper.sh --fast
```

## Package for preprints.org

Build a preprints-template preview PDF and source upload zip:

```bash
scripts/package_preprints.sh
```

Outputs:
- `paper/build/preprints_preview.pdf`
- `paper/build/preprints_source_upload.zip`
- `paper/archive/preprints_vN/current_source/`
- `paper/archive/preprints_vN/preprints_bundle/`

## Generate dashboard

```bash
TBD
```
