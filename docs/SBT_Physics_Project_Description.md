# Project Description

**A Physics is a Theory (Navier-Stokes Track)**  
SBT / "Six Birds" Follow-up  
Version 0.2 - 2026-02-11

## 1. Executive Summary
This repository is now focused on a single rigorous anchor: Navier-Stokes closure and blow-up analysis under the SBT framework. The objective is to study closure-consistent formulations, explicit assumption audits, and reproducible diagnostics that isolate where regularity arguments hold or fail.

## 2. Purpose and Positioning
The SBT math paper provides the abstract emergence calculus. This project applies that calculus to fluids with a proof-chain posture: formal theorem slots, explicit open lemmas, and reproducible numerical hooks that do not over-claim beyond proved assumptions.

## 3. Core Thesis
Navier-Stokes should be treated as a closure layer with explicit packaging, feasibility, and accounting constraints. Regularity claims must be conditional on auditable assumptions, and every hard step should be surfaced as a named settlement frontier.

## 4. Scope Boundaries and Non-goals
- No claim that the Clay problem is solved in this repository.
- No hidden substitution of numerics for proofs.
- No broad unification narrative beyond what is needed for the fluids theorem chain.

## 5. Primary Deliverables
- D1: Fluids theorem-chain documents (`docs/fluids/` and `docs/fluids/blowup/`).
- D2: Lean mechanization for algebraic/logic parts of the No-Zeno chain (`formal/Fluids/`).
- D3: Python diagnostics and stress sweeps for NS closure, route/capacity, and anti-localization (`python/src/nswave`, `python/scripts/ns*.py`, `python/scripts/nsbu*.py`).

## 6. Manuscript/Program Outline
- Closure-consistent NS setup and assumptions.
- Energy/accounting theorems and completion tracks.
- Blow-up chain scaffold with explicit hard lemma slots.
- Measurement hooks and falsification-oriented diagnostics.

## 7. Risks and Mitigations
| Risk | Failure mode | Mitigation |
| --- | --- | --- |
| R1 | Over-claiming from diagnostics | Keep all numeric evidence tagged as diagnostic-only unless proved. |
| R2 | Hidden BKM bottleneck reappears | Maintain explicit obstruction certificates and no-go statements in blow-up docs. |
| R3 | Assumption drift | Keep theorem inventory synchronized with code + Lean status. |

## 8. Active Documents
1. `docs/SBT_Physics_Step2_Fluids_NavierStokes_Rigorous.md`
2. `docs/fluids/theorem_inventory.md`
3. `docs/fluids/blowup/theorem_inventory_blowup.md`
4. `docs/fluids/blowup/clay_chain_scaffold.md`
