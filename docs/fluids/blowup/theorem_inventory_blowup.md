# Blow-up Track Inventory (NS)

## At-a-glance
- `MECH (Lean)`: abstract theorem-chain logic proved in Lean (`Zeno`, `ZenoWorkCap`, `RouteCapacity`, `ECTCap`, `AntiLocalization`, `IcapNoGo`).
- `NUMERIC (Python)`: simulation evidence and diagnostics only; supports plausibility, not theorem settlement.
- `OPEN (assumption)`: NS-specific hinge still unproved in this repository.

| ID | Claim (short) | Lean / Python anchor | Tag | Status |
|---|---|---|---|---|
| BU-T1 | Definitions + Zeno criterion + crossing-time logic | `formal/Fluids/Zeno.lean`; `python/scripts/nsbu01_dyadic_frontier_demo.py`; `python/scripts/nsbu02_latency_series.py` | `MECH (Lean)` | DONE |
| BU-T2 | Passivity + storage frontier => work quantum (`W^+ >= theta`) | `formal/Fluids/ZenoWorkCap.lean`; `docs/fluids/blowup/scale_passivity_axioms.md` | `MECH (Lean)` + `OPEN (assumption)` | PARTIAL |
| BU-T3 | Throughput + budget => latency bound (`Delta t >= theta/Cap`) | `formal/Fluids/ZenoWorkCap.lean`; `docs/fluids/blowup/obstruction_certificate.md` | `MECH (Lean)` + `OPEN (assumption)` | PARTIAL |
| BU-T4 | Divergence of `sum_j 1/Cap(j)` => No-Zeno | `formal/Fluids/Zeno.lean`; `formal/Fluids/ZenoWorkCap.lean` | `MECH (Lean)` | DONE |
| BU-T5 | ICAP-localization no-go pattern (discrete mechanized analogue) | `formal/Fluids/IcapNoGo.lean`; `docs/fluids/blowup/obstruction_certificate.md` | `MECH (Lean)` | DONE |
| BU-T6 | ECT finite-atom aggregation (`Lambda` grows with atom count) | `formal/Fluids/ECTCap.lean`; `formal/Fluids/RouteCapacity.lean` | `MECH (Lean)` | DONE (abstract) |
| BU-T7 | Conditional Clay chain assembly | `docs/fluids/blowup/clay_chain_scaffold.md` | `OPEN (assumption)` | PARTIAL |
| BU-H1 (HL-P2-ANTILOC) | Feasible inputs anti-localized (`eta_j -> 0`) | `formal/Fluids/AntiLocalization.lean` (discrete inequality); `python/scripts/nsbu28_antilocalization_diagnostic.py`; `python/scripts/nsbu29_antilocalization_sweep.py` | `MECH (Lean)` + `NUMERIC (Python)` + `OPEN (assumption)` | PARTIAL |
| BU-H2 (HL-P5-CHANNEL-DEL) | NS feasible set is finite-channel + delocalized (structural theorem) | `docs/fluids/blowup/antilocalization_hinge.md`; `python/scripts/ns_chan00_channel_count.py`; `python/scripts/ns_chan01_channel_delocalization.py` | `NUMERIC (Python)` + `OPEN (assumption)` | OPEN |
| BU-X1 | Toy Zeno gallery (necessity witness) | `python/scripts/nsbu12_toy_zeno_gallery.py` | `NUMERIC (Python)` | DONE |
| BU-HOOK-SGS | SGS bridge identification diagnostics | `python/scripts/nsbu04_sgs_frequency_response.py`; `python/scripts/nsbu05_sgs_realization_fit.py` | `NUMERIC (Python)` | HOOK |
| BU-HOOK-ROUTE | Route mismatch diagnostics | `python/scripts/nsbu06_route_dependence.py`; `docs/fluids/blowup/route_objects_ns.md` | `NUMERIC (Python)` | HOOK |

## Open gaps (NS-specific)
- Prove NS instantiation of ECT hypotheses (mode compression + per-atom ICAP uniformity).
- Prove HL-P2-ANTILOC structurally for NS feasible inputs, not only numerically.
- Prove HL-P5-CHANNEL-DEL (channel/addressability theorem for NS packaging).
- Prove HL-NS-ZENO bridge (finite-time blow-up implies Zeno frontier crossings).
- Resolve throughput/Cap assumptions without reintroducing the classical sup-norm obstruction channel.
