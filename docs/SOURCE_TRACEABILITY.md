# Numerical correction traceability

The mathematical generator is retained from `code_refined(8).m`. Original source archives and historical solver states remain unchanged in the separate private provenance deliverable. `STUDY_PROVENANCE.csv` retains original/release hashes for standalone sources and inputs. The archived 240-record reference CSV was not changed.

Changes from preparation commit `6f1b99847bc2f11580a7ef2c4ee1c20d14a42c7a`:

| File | Correction |
|---|---|
| `src/core/run_fixed_final_radius_VGPP.m` | Passive pre-canonicalization x/y capture and separate historical table output |
| `src/analysis/validate_historical_coordinates.m` | Verify geometry/statistical representations agree under endpoint mapping |
| `src/benchmark/run_benchmark_case.m` | Historical local-cell metrics; canonical alternative columns; separate geometry |
| `src/benchmark/benchmark_case_checkpoint.m` | Version 3, convention/metadata checks, preserve existing files |
| `config/vgpp_config.m` | Separate v3 checkpoint folder; numerical settings unchanged |
| `src/validation/migrate_v2_checkpoint.m` | Explicit restricted, lossless migration; no boundary-side inference |
| `src/validation/verify_reference_results.m` | Full-precision table support, maximum relative errors and mismatch counts; original absolute pass criteria |
| `src/validation/vgpp_release_checks.m` | Version/metadata/overwrite and historical-coordinate regressions |

Insertion, weighted legacy sampling, radius generation/growth, relaxation, seed assignments, stopping criteria and tolerances were not changed. Canonical particle tables matched saved configurations in each regenerated case checked by the running driver. Checkpoint source fingerprint is `8313894c98496ff31ef061f522dca73722ad27ff160105cedbec9ff86aa2c53b`. The fingerprint includes main/config/core/analysis/benchmark; validation utilities are outside this fingerprint and are covered by source inventory hashes.
