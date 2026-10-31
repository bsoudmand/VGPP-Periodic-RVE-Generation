# Numerical validation report — historical-binning correction

Preparation version `1.0.0-rc1`; status checked 2026-10-31T08:54:28.423154+00:00. **Classification B: Conditionally ready.** All 240 essential benchmark cases, the full coordinate/statistical recomputation, the standard output workflow and final MATLAB safeguard/dependency checks passed. Supplementary generator and continuation limitations remain explicitly identified; complete manuscript reproduction is not claimed.

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

## Formal numerical snapshot: 240/240 cases

Full-precision MATLAB checkpoint rows were compared with the matched archived records using `verify_reference_results`; all 240 identifiers are present and the full archived-reference gate passed. All 240 geometries, final radii and recomputed historical metrics passed an independent checkpoint audit. Four discrete fields match exactly; eleven continuous fields pass the original 1e-7 absolute tolerance. Relative discrepancies are descriptive, not replacement acceptance criteria; zero references use relative difference 0 for exact equality and infinity otherwise.

| Quantity | Max. absolute difference | Max. relative difference | Mismatches |
|---|---:|---:|---:|
| TargetN | 0 | 0 | 0 |
| ReachedTarget | 0 | 0 | 0 |
| NumberOfParticles | 0 | 0 | 0 |
| OverlappingPairs | 0 | 0 | 0 |
| FinalVF | 6.77236045021345e-15 | 1.22530823505845e-14 | 0 |
| DiscreteTargetVF | 4.9960036108132e-16 | 1.24892198405095e-15 | 0 |
| MeanRadius | 2.22044604925031e-16 | 1.2335811384724e-15 | 0 |
| RadiusScale | 2.22044604925031e-16 | 2.22044604925032e-16 | 0 |
| MaxOverlap | 0 | 0 | 0 |
| MinClearance | 4.43221848112074e-16 | 3.06561778277692e-13 | 0 |
| LocalVFRSD | 4.9960036108132e-16 | 4.43641728094663e-15 | 0 |
| LowVFCellFraction | 0 | 0 | 0 |
| MeanNND | 4.9960036108132e-16 | 1.38181880980074e-15 | 0 |
| NNDCV | 4.85722573273506e-16 | 4.78670977289517e-15 | 0 |
| RDFFirstPeak | 4.88498130835069e-15 | 2.99606701597357e-15 | 0 |


The compared dataset mixes compatible migrated records and regenerated cases. This is deterministic agreement to numerical tolerance for the checked seeded realizations, not merely agreement of ensemble distributions. All geometric and spatial descriptors were recomputed from the saved full-precision coordinates for all 240 audited cases, including NN/RDF measures, historical local statistics and canonical alternatives. Recomputed descriptors exactly matched the cached values (maximum difference 0). No bitwise equivalence to the original archived trajectories is claimed where original full trajectories are unavailable. Runtime is excluded from numerical equivalence testing.

See `BENCHMARK_SNAPSHOT.csv`, `CHECKPOINT_INTEGRITY.csv`, `BENCHMARK_PARTIAL_COMPARISON.csv`, and `BENCHMARK_CASE_STATUS.csv` (all 240 identities). A valid RSA-PBC realization below target is recorded separately from invalid geometry; censored supplementary runs are not geometric failures.

## Completed execution

At 2026-10-31T08:54:28.423154+00:00, 240 version-3 checkpoint files were saved; the serial corrected driver had visited 240/240 cases. The completed legacy producer had 240 saved cases and remains untouched. Its 240-case execution finished; the expected version-2 reference failure is confined to the two canonical local-bin descriptors. The run phase is `completed`. The full core-benchmark validation passed. Cases still lacking a version-3 checkpoint: 0. Exact identities appear in `BENCHMARK_CASE_STATUS.csv`.

Dense cases can require long execution; the earlier run includes a 33819.49 s measurement. All completed states and earlier incompatible-source evidence have been preserved. Source/configuration changes or another simultaneous generator launcher against the same v3 folder must be avoided. The final 240-case comparison, standard summaries and diagnostic figure-generation stages completed successfully.

## Remaining gates

All 240-case, full-reference, CSV/figure, and final targeted MATLAB checks passed; inspect the recorded scope before authorizing publication. The author-approved Figure 4 and its available R5-C13 source are preserved; its SD/truncation convention and numerical inputs are verified. Supplementary coverage and limits are in `MANUSCRIPT_STUDIES.md`; timing limits are in `RUNTIME_INVESTIGATION.md`. No complete manuscript reproduction is claimed. MIT/copyright approval, public visibility, GitHub release and public Zenodo deposit require owner authorization. No DOI, tag or release was created.

## Additional verification in this follow-up

The all-descriptor recomputation passed for 240 saved configurations, with source, complete configuration, format/convention, seed and continuous geometry checked. No saved audited geometry is invalid; no externally runtime-censored core case is recorded. Pending cases are not classified as failures or censored simulations. `BENCHMARK_CASE_STATUS.csv` identifies every intended case.

Cached-only execution of the original R1-C2 ordering script completed in MATLAB R2024b on all 18 expected coordinate caches, with generation and MC disabled in a scratch run-control copy. All 15783 scientific numerical values in the raw, summary, matched-effect, individual-RDF and RDF-summary tables agreed with the archive at 1e-7; maximum absolute difference was 4.61e-15. The two RDF files were byte-identical. Generation runtime and cache-load messages were excluded as execution metadata. This reproduces postprocessing, not the packing generators or MC trajectories.

The supplied approved Figure 4 PNG matches the current project export byte-for-byte. Its original 20-realization inputs and 72 exported statistics were independently checked; maximum absolute difference 5.12e-13. The manuscript defines SD bars and explains the triangle. This closes the earlier scientific/display-convention issue. Exact raster replay is not claimed.

Final audited execution origins: 136 losslessly migrated cases and 104 regenerated cases. Earlier version-2 states remain preserved.

The exported 240-row benchmark CSV was independently compared with the original unmodified reference CSV in Python. All 15 fields passed the same absolute thresholds; see `CSV_OUTPUT_COMPARISON.csv`. This supplements the full-precision MATLAB recomputation and the successful MATLAB output gate.

## Overlap reporting convention

The preserved legacy audit normalizes `MaxOverlap` to zero when the maximum overlap is no greater than 1e-8; `OverlappingPairs` counts tolerance violations. `MinClearance` retains the signed clearance. In the default zero-gap benchmark, the most negative clearance is -9.99990507066784e-9 (VGPP, target 0.65, run 8, seed 7808), so sub-tolerance overlap is present even though the reported `MaxOverlap` is zero. All 240 geometries meet the specified allowance; strict zero-overlap is not claimed. No output definition, reference value or tolerance was changed.
