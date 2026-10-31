# Numerical validation report — historical-binning correction

Preparation version `0.3.1-preview`; status checked 2026-10-31T07:03:59.196360+00:00. **Classification C: Not ready for public release.** The essential 240-case benchmark remains incomplete, and the approved manuscript Figure 4 error-bar provenance remains unresolved. These are completion/provenance gates; the corrected formal snapshot has no material reference mismatch.

## Executed verification

MATLAB 24.2.0.2712019 (R2024b), Windows 64-bit, with Image Processing Toolbox 24.2. Dependency analysis of `main.m` identifies MATLAB and Image Processing Toolbox. The targeted MATLAB regression suite passed: periodic minimum-image geometry; malformed/nonfinite input rejection; continuous overlap acceptance at 0.5 times and rejection at 1.5 times the 1e-8 tolerance; endpoint seed 7402; version-3 checkpoint reuse and identity; rejection of source/configuration/seed mismatch, version 2, missing historical metadata, and attempted overwrite with resume disabled; historical/canonical correspondence; legacy/log-key sampler numerical checks; optional LS analytical tests. Optional log-key selection remains a separate variant.

Both original R4-C10 and packaged R5-C13 archived-data figure scripts completed in MATLAB with identical exported statistics. Software OpenGL fallback was reported; figure generation completed. Earlier archived-data MATLAB geometry, 100-domain spatial, and LS consolidated audits remain applicable within their documented scope.

## Historical cell assignment

The original `code_refined(8).m` stores candidate-grid centers at either boundary, including exact L. The original local-cell function assigns x=L/y=L to the last cell when its half-open bin search finds no cell. Canonical x=0/y=0 goes to the first cell. Geometry is periodic-equivalent, but the cell statistics differ.

Production `run_fixed_final_radius_VGPP` now captures historical x/y immediately before endpoint canonicalization and returns a separate `historicalParticleTable`. This passive capture does not change radii, insertion, relaxation, stopping, or random draws. A source comparison confirmed that these capture/output statements are the only generator differences from the previously running corrected source. `run_benchmark_case` uses canonical coordinates for geometry, NN/RDF and area calculations, and historical coordinates for `LocalVFRSD` and `LowVFCellFraction`. Additional `CanonicalLocalVFRSD` and `CanonicalLowVFCellFraction` columns expose the representation-dependent alternative. `validate_historical_coordinates` checks identical particle identities/radii and exact coordinate correspondence under upper-endpoint mapping.

All 40 original low-density diagnostic cases were repeated through production, with seed and parameters unchanged. Historical local statistics passed the unchanged 1e-7 absolute tolerance, and canonical particle tables exactly matched the existing saved realizations. The original separate diagnostic had maximum historical RSD difference 4.44e-16 and zero low-fraction difference; production full-precision records are included in the formal comparison below. No reference values or tolerances were changed.

## Checkpoint audit and compatibility

Before integration, 204 version-2 cases were audited: every saved geometry, radius, source fingerprint and seed passed; no saved invalid geometry was found. Of these, 94 VGPP cases had canonical zero-boundary centers, making original boundary sides unavailable. At the subsequent migration scan, 206 saved cases were present: 111 were migrated losslessly and 95 required regeneration. These are different timestamped snapshots.

Version 3 stores historical and canonical tables, `historical-boundary-v1`, configuration and active-source SHA256. The default folder is `outputs/checkpoints-historical-v3`. Old files are neither overwritten nor deleted. An untouched 204-checkpoint pre-integration ZIP is retained locally. `migrate_v2_checkpoint` accepts only the explicitly audited prior fingerprint and identical numerical configuration; VGPP cases with any zero-boundary center are rejected for migration and regenerated. No zero-to-L inference is performed. Eligible saved cases retain their measured generation times, with migration metadata identifying their origin. The migration is explicit, not automatic acceptance by the resume guard.

Old fingerprint: `80af3a550fdec5744dbfbb3934fe3ed6db2fc532edbfa6f4591f7ef5931cc391`.

Current active fingerprint: `8313894c98496ff31ef061f522dca73722ad27ff160105cedbec9ff86aa2c53b`.

## Formal numerical snapshot: 183/240 cases

Full-precision MATLAB checkpoint rows were compared with the matched archived records using `verify_reference_results`; this partial check uses only the matched reference subset and does not satisfy the full 240-record gate. All 183 geometries, final radii and recomputed historical metrics passed an independent checkpoint audit. Four discrete fields match exactly; eleven continuous fields pass the original 1e-7 absolute tolerance. Relative discrepancies are descriptive, not replacement acceptance criteria; zero references use relative difference 0 for exact equality and infinity otherwise.

| Quantity | Max. absolute difference | Max. relative difference | Mismatches |
|---|---:|---:|---:|
| TargetN | 0 | 0 | 0 |
| ReachedTarget | 0 | 0 | 0 |
| NumberOfParticles | 0 | 0 | 0 |
| OverlappingPairs | 0 | 0 | 0 |
| FinalVF | 6.55031584528842e-15 | 1.20510836638324e-14 | 0 |
| DiscreteTargetVF | 4.9960036108132e-16 | 1.24892198405095e-15 | 0 |
| MeanRadius | 2.22044604925031e-16 | 1.2335811384724e-15 | 0 |
| RadiusScale | 2.22044604925031e-16 | 2.22044604925032e-16 | 0 |
| MaxOverlap | 0 | 0 | 0 |
| MinClearance | 4.43221848112074e-16 | 3.06561778277692e-13 | 0 |
| LocalVFRSD | 4.9960036108132e-16 | 4.43641728094663e-15 | 0 |
| LowVFCellFraction | 0 | 0 | 0 |
| MeanNND | 4.9960036108132e-16 | 1.34639957613479e-15 | 0 |
| NNDCV | 4.85722573273506e-16 | 4.78670977289517e-15 | 0 |
| RDFFirstPeak | 4.88498130835069e-15 | 2.99606701597357e-15 | 0 |


The compared dataset mixes compatible migrated records and regenerated cases. This is deterministic agreement to numerical tolerance for the checked seeded realizations, not merely agreement of ensemble distributions. Other spatial descriptors in migrated records retain their audited prior generation outputs; geometry and historical local descriptors were recomputed during the new checkpoint audit. No bitwise equivalence to the original archived trajectories is claimed where original full trajectories are unavailable. Runtime is excluded from numerical equivalence testing.

See `BENCHMARK_SNAPSHOT.csv`, `CHECKPOINT_INTEGRITY.csv`, `BENCHMARK_PARTIAL_COMPARISON.csv`, and `BENCHMARK_CASE_STATUS.csv` (all 240 identities). A valid RSA-PBC realization below target is recorded separately from invalid geometry; censored supplementary runs are not geometric failures.

## Live progress, separate from the formal snapshot

At 2026-10-31T07:03:59.196360+00:00, 184 version-3 checkpoint files were saved; the serial corrected driver had visited 156/240 cases. Migrated future-target cases explain why these counts differ. The legacy producer had 212 saved cases and remains untouched. The run phase is `running`. No full-validation pass is asserted. Cases still lacking a version-3 checkpoint: 56. Exact identities appear in `BENCHMARK_CASE_STATUS.csv`.

Dense cases can require long execution; the earlier run includes a 33819.49 s measurement. Ongoing jobs and checkpointed progress have been preserved rather than declaring completion. Source/configuration changes or another simultaneous launcher against the same v3 folder must be avoided. Final 240-case comparison, standard summary generation and new-run figures remain pending.

## Remaining gates

Complete and audit all 240 cases; inspect all completion statuses; execute the final full-reference and CSV/figure gates. Resolve or explicitly identify the missing manuscript Figure 4 error-bar source. Supplementary coverage and limits are in `MANUSCRIPT_STUDIES.md`; timing limits are in `RUNTIME_INVESTIGATION.md`. No complete manuscript reproduction is claimed. MIT/copyright approval, public visibility, GitHub release and public Zenodo deposit require owner authorization. No DOI, tag or release was created.
