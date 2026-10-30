# Technical validation report — release preparation

Repository: https://github.com/bsoudmand/VGPP-Periodic-RVE-Generation  
Preparation version: `0.3.1-preview` (not a published release).  
Compared checkpoint snapshot: 2026-10-30T20:20:35+00:00.

## Release decision

**NOT READY for v1.0.0 or a permanent scientific release.** The full 240-case run has been started but is incomplete, and the available rerun records contain two classes of local-bin statistic discrepancies. Public visibility, MIT adoption, GitHub release publication, and Zenodo publication remain unapproved. No tag or DOI is assigned by this preparation.

## Benchmark execution and source changes

MATLAB `24.2.0.2712019 (R2024b)` with Image Processing Toolbox executed the tests. Dependency analysis with the active folders on the MATLAB path identified MATLAB and Image Processing Toolbox. The run uses Lx=Ly=10, radius 0.18, grid 250-by-250, overlap tolerance 1e-8, the six original target fractions, 20 realizations per method/fraction, the original seed policy, the legacy weighted sampler, and unchanged convergence/stopping settings.

The first full attempt on the previous modular code stopped at VGPP seed 7402 after saving the first VGPP/RSA pair. The inclusive candidate grid could return a center exactly at L, which the newly strict canonical-coordinate audit rejected. Those original-attempt checkpoints and logs were retained locally.

The candidate maps only coordinates exactly equal to Lx/Ly to zero immediately before its final audit. It does not map arbitrary out-of-domain or nonfinite values. Sampling, insertion, growth, pair distances, and tolerances were not changed. A separate checkpoint correction requires the stored seed to match `7300 + 100*targetIndex + runIndex`. New source identifiers prevent reuse of checkpoints made with the old code. A reusable targeted regression function was added under `src/validation/`.

Active benchmark source fingerprint: `80af3a550fdec5744dbfbb3934fe3ed6db2fc532edbfa6f4591f7ef5931cc391`. The corrected run started from fresh checkpoints; the incompatible original-attempt pair was not silently reused. A duplicate launcher briefly overlapped the corrected run and was stopped. One completed checkpoint identity is retained per method/fraction/run/seed, and the snapshot contains no duplicate identities. New runtimes also reflect concurrent audit work; they are not isolated performance measurements.

## Executed checks

| Check | Result | Scope |
|---|---|---|
| Core `vgpp_self_test` | PASS | Periodic geometry, malformed inputs, target-count convention, sampler tests |
| Continuous-coordinate 1e-8 tolerance boundary | PASS | Permitted overlap at 0.5 times tolerance; rejection at 1.5 times tolerance |
| Invalid coordinates/radii and array shapes | PASS | NaN/Inf, nonpositive radii, noncanonical centers, shape mismatch |
| Seed 7402 endpoint regression | PASS | 393-particle target attained with canonical coordinates |
| Checkpoint reuse without duplicate case/write | PASS | Same record returned, same checkpoint size/modification time, one case file |
| Checkpoint parameter/source/seed mismatch rejection | PASS | Incompatible metadata rejected with checkpoint-mismatch error |
| Optional LS analytical tests | PASS | Collision, periodic-image, impulse, audit, and heap tests; not full LS validation |
| Original full-precision coordinate audit | PASS | Saved 747-particle high-density configuration, allowance 1e-8 |
| Original independent domain-size spatial audit | PASS | 100 archived coordinate configurations; no generator rerun |
| Original consolidated LS archive audit | PASS | Saved production results/states audited read-only on scratch copies |
| Archived Figure 4 postprocessing | PASS | PNG and calculated CSV generated from the original archived 240 records; not regenerated simulation results |
| Full 240-case end-to-end reference validation | INCOMPLETE / discrepancies found | No full pass is claimed; new-run figure generation remains gated |

The source and parameter mismatch tests reject changed identifiers; they do not constitute a guarantee against arbitrary manual checkpoint corruption. The main checkpoint format stores complete parameters, a format version, source SHA-256, case metadata, and particle coordinates.

## Quantitative partial comparison

136 unique checkpoint cases were inspected, with passing coordinate/radius audits and matching source/seed identifiers. Four discrete fields were compared exactly and eleven continuous fields at absolute tolerance 1e-7. Wall-clock time was excluded from numerical equivalence tests. See `BENCHMARK_SNAPSHOT.csv`, `BENCHMARK_PARTIAL_COMPARISON.csv`, and `CHECKPOINT_INTEGRITY.csv`.

| Field | Compared cases | Mismatches | Maximum absolute difference |
|---|---:|---:|---:|
| TargetN | 136 | 0 | 0 |
| ReachedTarget | 136 | 0 | 0 |
| NumberOfParticles | 136 | 0 | 0 |
| OverlappingPairs | 136 | 0 | 0 |
| FinalVF | 136 | 0 | 7.1054273576e-15 |
| DiscreteTargetVF | 136 | 0 | 0 |
| MeanRadius | 136 | 0 | 0 |
| RadiusScale | 136 | 0 | 0 |
| MaxOverlap | 136 | 0 | 0 |
| MinClearance | 136 | 0 | 4.39969241595e-16 |
| LocalVFRSD | 136 | 34 | 0.0140414990048 |
| LowVFCellFraction | 136 | 13 | 0.015625 |
| MeanNND | 136 | 0 | 0 |
| NNDCV | 136 | 0 | 1.11022302463e-16 |
| RDFFirstPeak | 136 | 0 | 0 |

The two failing fields are center-based local cell-bin statistics. For seed 7402, local RSD changes from 0.196509342287901 to 0.198639948000739 while particle count, area fraction, clearance, NN-CV, and RDF peak agree. Moving a periodic endpoint from L to 0 leaves the torus geometry unchanged but moves the center from the legacy last bin to the first bin. The reference therefore cannot be declared reproduced under the canonical storage convention. No sampler switch, relaxed tolerance, or modified reference value was used to conceal this discrepancy. Its scientific reconciliation remains a release blocker.

Partial timing distributions for matched cases are recorded in `BENCHMARK_RUNTIME_PARTIAL.csv`. Compare attempted and successful populations consistently. Different hardware and concurrent workload affect elapsed time; exact equality is not expected, and performance equivalence has not been established.

## Supplementary evidence

All ten requested study families have original source implementations and archived outputs available in the project; 33 selected standalone scripts are organized in the repository. The clean supplementary archive contains 268 CSV/settings entries. Original source hashes and I/O-only path normalization are recorded in `STUDY_PROVENANCE.csv`. No numerical kernel was merged or replaced.

Independent pairwise checks covered 182 archived coordinate files: none exceeded its applicable 1e-8 overlap/gap allowance; all numerical inputs/radii were finite and positive. Exact upper-endpoint coordinates occur in 100 files. Twelve files do not have all radii equal to 0.18, including censored growing-disk states and a separate non-baseline configuration; they must not be treated as attained fixed-final-radius benchmark cases. The 0.76 configuration contains 985 strictly overlapping pairs, but its maximum overlap, 9.9982804813159e-9, remains within the stated numerical allowance. Strict nonintersection is not claimed.

The optimized LS production dataset records 7 completed cases and 2 event-budget-censored cases out of 9; only 1 of 3 at 0.74 completed. The older LS dataset contains 9 wall-budget-censored cases. The focused ablation data contain 12 successes, 3 completed-below-target cases, and 3 runtime-censored cases. Censoring is not a packing failure. See `MANUSCRIPT_STUDIES.md` for exact script/data mappings, execution instructions, missing per-snapshot evidence, and continuation-state requirements.

## Remaining work and publication gates

- Finish or safely resume the remaining benchmark cases; retain checkpoints and examine the complete record/geometry comparison.
- Resolve and document the legacy-versus-canonical local-bin statistic discrepancy before claiming reference equivalence or preparing v1.0.0. This is not a hardware-timing difference.
- Complete any manuscript-critical supplementary verification not covered by archived-data audits, including ordering/MC postprocessing and required exact continuation/figure inputs.
- Verify the full main/plot pipeline after a scientifically justified reference decision; archived Figure 4 execution alone is insufficient.
- Obtain explicit owner approval for the MIT terms/copyright line and public visibility; obtain separate authorization for release/Zenodo publication.

No complete manuscript reproducibility, unbiased legacy sampling, equilibrium MC ensemble, LS packing ceiling, or permanent archive is claimed.
