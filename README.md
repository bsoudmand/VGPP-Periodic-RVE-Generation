# Void-Guided Particle Packing for Periodic RVE Generation

MATLAB implementation of Void-Guided Particle Packing (VGPP) for generating two-dimensional periodic particulate representative volume elements (RVEs) composed of circular particles. The consolidated workflow implements a fixed-final-radius, monodisperse benchmark and comparison with random sequential adsorption under periodic boundary conditions (RSA-PBC).

**Version:** `0.3.1-preview`  
**Repository:** https://github.com/bsoudmand/VGPP-Periodic-RVE-Generation

## Research purpose and validation status

The software supports investigation of attainable particle area fraction, geometric admissibility, computational cost, and spatial statistics in periodic disk packings. These geometric configurations are intended for particulate microstructure studies; the code does not demonstrate that an individual generated domain is statistically representative for a particular material property.

**Full validation has not passed.** A checkpointed 240-case rerun is in progress. At the numerically inspected snapshot (2026-10-30T20:33:45+00:00), 162 unique cases were available. Counts, target-attainment flags, radii, periodic-overlap measures, NN statistics, and RDF first peaks agreed with the matched reference records within the configured tolerances. However, `LocalVFRSD` and `LowVFCellFraction` disagreed in some cases because exact upper-boundary centers are now stored at their equivalent zero-boundary positions. The legacy cell-binning calculation is sensitive to this representation. The archived reference CSV and numerical tolerances have not been altered to force agreement.

MATLAB R2024b core/safeguard regressions, checkpoint compatibility/resume tests, and optional LS analytical tests passed. Three original archived-data MATLAB audits also completed. These checks do not establish full manuscript reproducibility or authorize a public release. See `docs/VALIDATION_REPORT.md` for quantitative results and `docs/MANUSCRIPT_STUDIES.md` for study coverage and limitations.

## Main algorithm features

- Periodic occupancy and distance-transform maps for void-guided candidate generation, including largest-void and weighted candidate selection.
- Batch insertion with periodic minimum-image geometric checks.
- Temporary reduction of seed radii followed, where required, by gradual growth and collective repulsive rearrangement toward the prescribed final radius.
- Stochastic relaxation during seed construction and a final periodic pairwise overlap audit.
- A reduced-count recovery procedure when the requested fixed-radius packing cannot be attained; target attainment and geometric validity are recorded separately.
- RSA-PBC comparison, nearest-neighbor and radial-distribution statistics, and local area-fraction measures.
- Per-realization checkpoints with configuration and active-source fingerprint checks.

The default workflow uses no ordered or hexagonal fallback. A requested target density is not guaranteed to be attained. Geometrically invalid cases are recorded as failures and their spatial statistics are excluded from valid-result summaries.

## Software requirements

- MATLAB **R2024b**; compatibility with other releases has not been established for this package.
- **Image Processing Toolbox**, required for `bwdist` and `imregionalmax`.
- MATLAB's standard table, plotting, and file I/O functionality, and a working Java runtime for the active-source SHA-256 fingerprint.
- Write permission in the local project/output directories.

The default workflow does not require Parallel Computing Toolbox. MATLAB and its required toolbox must be obtained separately; they are not distributed with this repository.

## Repository structure

```text
main.m                         Benchmark entry point
config/vgpp_config.m           Configuration and output paths
src/core/                     Packing, periodic geometry, relaxation, audit
src/benchmark/                Drivers, checkpoints, figures, source fingerprint
src/analysis/                 Spatial statistics and reference summaries
src/validation/               Self-tests and reference-comparison routine
inputs/reference/             Archived 240-record benchmark CSV
studies/optional/             Modular LS reference and analytical tests
studies/manuscript/           Original standalone manuscript study families
inputs/supplementary_reference.zip  Archived study CSVs, coordinates, settings
outputs/benchmark/            Generated numerical tables
outputs/figures/              Generated figures
outputs/checkpoints/          Generated per-realization MAT checkpoints
docs/                         Provenance, kernel audit, scope, release checklist
CITATION.cff                  Software citation metadata
LICENSE_SELECTION_REQUIRED.md License decision pending
```

The checkpoint directory is created during execution. Generated outputs are ignored by Git, except for the directory placeholders.

## Installation

Clone the repository, or download its ZIP archive and extract it to a writable directory. While the repository is private, access requires authorization from its owner.

```bash
git clone https://github.com/bsoudmand/VGPP-Periodic-RVE-Generation.git
```

Open MATLAB R2024b, confirm that Image Processing Toolbox is installed using `ver`, and select the extracted or cloned project root as the current folder. `main.m` adds the active source folders to the MATLAB path.

## Running the workflow

From the project root, run:

```matlab
main
```

Alternatively, open `main.m` in the MATLAB Editor and select **Run**. The script clears the workspace and closes existing figures before configuring the benchmark. It runs the included self-tests, then executes or resumes the target-area-fraction sweep.

The default study comprises six nominal area fractions, 20 realizations per fraction, and two methods: 240 attempted method/fraction/realization cases. Dense-packing cases can be computationally expensive. No runtime guarantee is provided.

For an initial smoke test, use a separate project copy and change the following setting in `config/vgpp_config.m`:

```matlab
params.targetVFSweepRuns = 2;
```

Then run `main`. A two-realization smoke test is exploratory and does not run the full archived-reference comparison. The self-test contains assertions tied to the default geometry, radius, and tolerance; changing these settings may require corresponding independent test changes and validation.

## Input parameters

Edit `config/vgpp_config.m` before execution. All lengths use a consistent, user-defined model-length unit; the configuration assigns no physical unit.

| Parameter | Default | Role |
|---|---|---|
| `Lx`, `Ly` | `10`, `10` | Periodic domain dimensions |
| `Nx`, `Ny` | `250`, `250` | Occupancy/distance-transform grid dimensions |
| `rMono` | `0.18` | Prescribed final disk radius |
| `particleMode` | `"monodisperse"` | Default particle-size setting |
| `minGap` | `0` | Nominal minimum surface clearance |
| `targetVFSweepList` | `[0.40 0.50 0.55 0.60 0.65 0.70]` | Nominal particle area fractions |
| `targetVFSweepRuns` | `20` | Realizations per target and method |
| `runRSAInTargetSweep` | `true` | Include RSA-PBC comparison |
| `batchSize`, `nCandidates` | `80`, `4000` | Insertion batch and candidate counts |
| `distanceUpdateInterval` | `5` | Distance-map refresh interval |
| `useLargestVoidInsertion`, `largestVoidFraction` | `true`, `0.70` | Largest-void candidate selection |
| `weightedSamplingAlgorithm` | `'legacy'` | Default sampling rule retained from the source benchmark |
| `fixedRadiusSeedVF` | `0.50` | Seed-packing radius-reduction control |
| `fixedRadiusGrowthScale` | `1.0025` | Trial multiplicative radius-growth factor |
| `fixedRadiusMaxGrowthSteps` | `300` | Growth-step limit |
| `maxCycles`, `maxFailedCycles` | `10000`, `500` | Seed insertion stopping limits |
| `rsaMaxTrials`, `rsaMaxFailedTrials` | `5000000`, `200000` | RSA-PBC stopping limits |
| `enableCheckpoints`, `resumeCheckpoints` | `true`, `true` | Save and reuse compatible completed cases |

Additional relaxation, backtracking, candidate-separation, spatial-analysis, tolerance, and path settings are documented by their fields in the configuration. The presence of radius-distribution fields does not establish a validated polydisperse benchmark: `main.m` uses the fixed-final-radius workflow.

The particle count is rounded to an integer:

```matlab
Ntarget = max(1, round(phiA * Lx * Ly / (pi * rMono^2)));
```

Consequently, the discrete target area fraction can differ from the nominal value. In this two-dimensional model, particle fraction denotes **area fraction**, including where legacy field or function names contain `VF` or `volume_fraction`.

## Expected outputs

Successful completion of the relevant stages produces:

- `outputs/benchmark/Table_R1_fixedRadius_raw.csv`: per-case counts, seeds, target attainment, geometry status, runtime, and spatial statistics.
- `outputs/benchmark/Table_R2_fixedRadius_summary.csv`: method/target summary statistics.
- `outputs/checkpoints/case_*.mat`: saved case records, particle tables (`ID`, `x`, `y`, `radius`), configuration, and source fingerprint.
- `outputs/benchmark/reference_comparison.csv`: full-benchmark comparison report, if the comparison stage is reached.
- `outputs/figures/`: benchmark figures supported by the included plotting functions, subject to successful completion of preceding stages.
- For the default 20-realization comparison workflow, `Fig4_ordering_runtime_original20.png` and `Fig4_original20_values.csv` are generated after successful reference validation. The runtime plot uses the documented R5-C13 clipped-whisker convention; consult the plotting source and exported statistics when interpreting it.

The supplied `inputs/reference/Table_R1_fixedRadius_raw.csv` is an archived input and is not overwritten. Completion states distinguish `target_attained`, `valid_below_target`, and `invalid_geometry`. An invalid geometry is not a successful RVE; affected area-fraction/spatial quantities are stored as `NaN` where implemented. Aggregate area-fraction and spatial summaries must be interpreted together with attainment and validity status.

## Numerical tolerances

| Setting or criterion | Value | Interpretation |
|---|---|---|
| `overlapAuditTolerance` | `1e-8` | Permitted pairwise overlap in model-length units in the geometric audit |
| `growthCompressionOverlapTolerance` | `1e-8` | Configured overlap tolerance for growth/relaxation; coupled to the audit tolerance |
| `targetVFTolerance` | `1e-6` | Configured area-fraction tolerance; the fixed-radius driver determines attainment using final-radius validity and the integer target count |
| `validationTolerance` | `1e-7` | Absolute tolerance for finite, continuous metrics in the archived-reference comparison |
| Final radius consistency | `1e-12` | Absolute radius check used in the VGPP geometry-validity flag |

Reference identifiers and discrete counts/status fields are compared exactly. Nonfinite-value patterns must also match. Elapsed runtime is excluded from numerical equivalence testing. The overlap allowance does not imply strictly disjoint disks in exact arithmetic, and the distance-transform grid is a candidate-generation discretization rather than the final overlap acceptance criterion.

## Reproducibility and validation procedure

1. Record the Git commit, MATLAB release, toolbox version, operating system, hardware, and complete configuration. Retain an unmodified copy of the archived input CSV.
2. Run the two-realization smoke test in a separate copy. Inspect geometry/failure records and checkpoint creation. Rerun the same configuration to assess resume behavior. Do not treat this as a full numerical validation.
3. For the reference study, restore the default settings and use a new, empty checkpoint/output location. Use 20 realizations at the six supplied area fractions with RSA-PBC enabled.
4. The per-case seed is `7300 + 100*targetIndex + runIndex`, reset separately before each method. Preserve target ordering, run numbering, and the default `'legacy'` sampling rule when investigating agreement with the archived study.
5. Inspect every validity and completion status. `verify_reference_results` aligns records by method, target fraction, run, and seed; checks identifiers/counts, finite metrics, and nonfinite patterns; and rejects rerun records flagged as geometrically invalid.
6. In the default 20-realization/RSA workflow, `main.m` stops before figure generation if the reference is missing or the comparison fails. Preserve the diagnostic report when produced and investigate all discrepancies. A failed or incomplete comparison is not validated reproduction.
7. Checkpoints require matching configuration, format version, and active-source SHA-256 fingerprint. Source/configuration changes require separate outputs; do not mix checkpoints from different experiments.
8. The opt-in `'log_keys'` sampler addresses a numerical-stability issue in the legacy power-based selection. It can alter candidate selection and packing outcomes and requires a separate benchmark. The included sampler self-test does not establish distributional equivalence or unbiasedness of the legacy sampler.

## Scope and known limitations

The automated core benchmark covers two-dimensional, monodisperse, fixed-final-radius disks under periodic boundary conditions. Original supplementary study families and archived numerical reference files are additionally provided under `studies/manuscript/` and `inputs/supplementary_reference.zip`; they are not invoked by `main.m`. Their kernels and sampler choices can differ from the core benchmark and have not been silently unified.

The study inventory covers high-density packing, matched-physics comparator/ablation studies, LS, grid/gap/domain sensitivity, ordering and finite-chain MC, and runtime/scaling. Independent geometry checks covered 182 archived coordinate files; full simulation reruns and all figure workflows have not been verified. Some exact LS continuation and MAT re-export workflows require original solver states preserved in the separate provenance deliverable, not in the clean release. MC snapshot metrics are present but full coordinates for every snapshot are not.

Two limited code corrections are documented: exact upper-endpoint output coordinates are mapped to zero without wrapping arbitrary invalid coordinates, and checkpoint reuse additionally verifies the stored seed. The legacy weighted sampler, insertion/growth rules, convergence settings, and tolerances remain unchanged. The boundary correction changes legacy local cell-bin statistics in some low-density records, so full agreement with the original archive is not claimed.

To run the targeted regression suite after adding `src/validation` to the path:

```matlab
vgpp_release_checks(pwd)
```

This suite exercises geometry/tolerance safeguards, the seed-7402 endpoint regression, compatible checkpoint reuse, rejection of mismatched metadata, and optional LS analytical tests. It does not replace the 240-case comparison. See `docs/MANUSCRIPT_STUDIES.md` for standalone study installation, data staging, required inputs, and execution instructions. A `v1.0.0` release candidate must wait until the unresolved comparison discrepancies and remaining release checks have been resolved.

## Citation and license

Use `CITATION.cff` to cite the software authors, title, version, repository, and the exact commit used. Cite any associated article separately once its bibliographic details are verified. No publication or software DOI is assigned in this metadata.

An open-source license has not yet been approved. MIT is proposed, but no MIT permission grant is applied by this README. The repository remains private pending the owner's approval for public visibility. See `LICENSE_SELECTION_REQUIRED.md`.
