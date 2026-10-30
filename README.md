# Void-Guided Particle Packing (VGPP)

MATLAB R2024b implementation of the **fixed-final-radius, 2D periodic disk** benchmark. The entry point is `main.m`. It generates VGPP and RSA-PBC realizations at nominal particle area fractions 0.40, 0.50, 0.55, 0.60, 0.65, and 0.70, using the same prescribed final particle radius and periodic minimum-image geometric acceptance criterion.

## Requirements

- MATLAB **R2024b**.
- **Image Processing Toolbox** for `bwdist` and `imregionalmax`.
- Standard MATLAB functionality for tables, file export, and plots.
- The optional LS reference (`studies/optional`) uses standard MATLAB functions, but is a research benchmark with separate checkpoints and validation needs.

## Execution

1. Unzip the repository to a folder with write permission.
2. From MATLAB, open the folder and run `main` (or `run('main.m')`).
3. The default **full benchmark** uses 20 fixed seeds per target per method. Some dense-packing cases are computationally expensive and can take many hours. Use `params.targetVFSweepRuns = 2` in `config/vgpp_config.m` only for a preliminary smoke test; do not compare two-run results against the full reference.
4. The benchmark writes CSV files under `outputs/benchmark/` and figures under `outputs/figures/`. Existing results in `inputs/reference/` are never overwritten.

The seed policy is retained from the original benchmark: `seedNow = 7300 + 100*targetIndex + runIndex`, resetting `rng(seedNow)` before each VGPP and RSA-PBC realization. The target particle count is `round(phiA*Lx*Ly/(pi*rFinal^2))`.

## Directory layout

| Directory | Purpose |
|---|---|
| `main.m` | One-command benchmark entry point |
| `config/` | Fixed benchmark settings |
| `src/core/` | Independent generation, PBC geometry, relaxation and audit functions |
| `src/benchmark/` | Original multi-target driver and figure-export helpers |
| `src/analysis/` | Spatial descriptors and archived-data summaries |
| `src/validation/` | Geometric smoke tests and record-by-record comparisons |
| `inputs/reference/` | Supplied original benchmark CSV, used for comparisons |
| `studies/optional/` | Optional independent LS reference (not part of default workflow) |
| `outputs/` | Generated files; safe to delete and rerun |
| `docs/` | Source inventory, unresolved dependencies and audit notes |

## Scientific safeguards

- Final particle radius `r=0.18`; reduced seed radii are enlarged only up to this value.
- Minimum-image distances are used for periodic acceptance.
- Final overlap acceptance tolerance: `1e-8` model-length units, independent of distance-transform grid resolution.
- No ordered/hexagonal fallback is enabled in the benchmark.
- A failed target is recorded separately from a geometrically accepted reduced-count configuration. Geometrically invalid configurations are flagged and excluded from spatial result summaries.
- RSA-PBC spatial comparisons should include only attained target densities; the reachability plot includes all attempted targets.
- Elapsed time depends on hardware and workload; it is **not** included in exact record-by-record numerical validation.
- The original benchmark plots are retained. The full 20-run workflow additionally produces `Fig4_ordering_runtime_original20.png` using the R5-C13 clipped-whisker convention and exports the actual statistics as `Fig4_original20_values.csv`. Other later layout-only manuscript panels are not reconstructed by this package.

## Resume and failure classification

The benchmark saves each method/target/run immediately to `outputs/checkpoints/` using an atomic temporary-file rename. Rerunning `main` resumes compatible checkpoints. Each MAT file carries a format version and a SHA-256 fingerprint of the active `main.m`, `config/`, and MATLAB functions in `src/core/`, `src/benchmark/`, and `src/analysis/`. Any change to those source bytes or the numerical configuration prevents a checkpoint from being reused, even when the parameter values are unchanged. The fingerprint is computed using MATLAB’s Java runtime. To start an independent experiment, use a separate repository copy or explicitly clear the checkpoint folder. Full results are assembled only after every requested case has completed.

Raw results include `GeometryValid` and `CompletionStatus` (`target_attained`, `valid_below_target`, `invalid_geometry`). Invalid configurations are retained only as traceable audit records; their spatial metrics and achieved area fraction are stored as `NaN`. These classification changes do not deliberately alter particle-generation mathematics.

## Weighted sampler modes

`params.weightedSamplingAlgorithm = 'legacy'` preserves the original published numerical calculation, including potential underflow in its power-based keys. Set it to `'log_keys'` for an **experimental numerical-stability correction**. That change can alter candidate selection and the eventual packing output; it must be benchmarked and validated independently before adoption in the published reference workflow. The self-test draws 4,000 out of 50,000 equal-weight candidates with both algorithms, checks uniqueness, and checks that the log-key selection covers the index range. It does not prove numerical equivalence or statistical unbiasedness of the legacy method.

## MATLAB validation protocol

1. Run `vgpp_self_test(vgpp_config(pwd))` after adding `config/`, `src/core/`, and `src/validation/` to the path, or simply start with `main` (which invokes the self-test).
2. Set `targetVFSweepRuns=2` for a smoke test in a **separate copy** of the repository. Confirm checkpoint creation, and rerun `main` to verify successful resume.
3. For the original seed policy, configure `targetVFSweepRuns=20`; use an empty checkpoint directory. Check the console for any `invalid_geometry` cases, inspect raw `GeometryValid`, and compare exact target counts and geometry metrics with the archived file.
4. Review `reference_comparison.csv`, including the `NonfinitePatternMatch` and `Pass` columns. The full numerical comparison must pass before final manuscript-figure generation. `main.m` stops before figure output if this check fails. Elapsed runtimes depend on hardware.
5. Evaluate `'log_keys'` only in a separate clean run. Do **not** mix these checkpoints with `'legacy'`.

## Numerical safeguards added after the second review

- The final geometric audit explicitly rejects nonfinite center coordinates, nonpositive or nonfinite radii, unequal particle-array lengths, and particle centers outside the canonical periodic interval `[0,L)`. Noncanonical coordinates raise a MATLAB error rather than silently passing an overlap comparison.
- Existing checkpoints from earlier previews have a different format and must not be resumed. Start with an empty `outputs/checkpoints/` directory or use a separate repository copy.
- For the complete 20-seed VGPP/RSA-PBC benchmark, numerical reference checking now precedes all publication-figure generation. A failed comparison leaves its diagnostic CSV and aborts plotting. Smoke-test figures (2-run setting) are exploratory and are not reference-validated.
- The active source checksum covers executable benchmark code and configuration, not the optional LS study or reviewer-specific scripts outside this package.

## License and scope

`LICENSE_SELECTION_REQUIRED.md` is an explicit release blocker. No permission-granting LICENSE is included because the authors have not selected one. Choose an appropriate license (e.g., MIT or BSD-3-Clause), approve its terms, and add its full text as `LICENSE` before uploading publicly. The package is the fixed-final-radius VGPP/RSA-PBC benchmark only; reviewer-specific investigations are not yet fully integrated.

## Release checklist

Review `docs/RELEASE_CHECKLIST.md` before public release. The reviewer-specific experiments are still outside the automated main workflow, and MATLAB execution remains outstanding.

## Validation status

The code has been **statically inspected and refactored**, but MATLAB R2024b was **not available** in the creation environment. Therefore, the modular version has **not yet passed an actual MATLAB execution or numerical rerun**. `verify_reference_results` compares rerun CSVs to the supplied original 240-row benchmark, records differences, and emits an explicit warning if they exceed tolerances. A static code review does not substitute for execution tests.

Further reviewer-specific experiments (high-density ceiling, varied-clearance trials, bidisperse cases, LS continuation, finite-duration MC, and large robustness/scaling studies) are **not equivalent to a single benchmark rerun**. Several use distinct kernels or depend on preserved coordinates/checkpoints. Their original versions are in the separately supplied provenance archive and must be integrated only after a matched validation run. Do not treat this core package as a complete numerical reproduction of every supplementary result.
