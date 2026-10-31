> Historical preparation snapshot. The current correction and executed results are documented in `VALIDATION_REPORT.md` and `SOURCE_TRACEABILITY.md`.

# Response to independent repository review

This document records changes made after the second static code review. None of the modifications below constitutes MATLAB execution validation.

| Finding | Action | Status |
|---|---|---|
| Power-key underflow in weighted sampling | Added opt-in `log_keys` implementation using `log(U)/weight`; retained original `legacy` mode as default. Threaded configuration through both sampling selectors. | Both modes require MATLAB execution; numerical equivalence is **not** asserted. |
| No crash recovery | Added a per-method, per-target, per-seed MAT checkpoint, saved through temporary-file replacement. On resume, the complete configuration and run identifiers must match. | Static implementation only; resume behavior must be tested in MATLAB. |
| Geometrically invalid outputs | Added `GeometryValid`, `CompletionStatus`; invalid geometry gets `NaN` in spatial outcome statistics and attained fraction; audits and elapsed time remain in raw records. | Geometry must be checked with MATLAB tests and benchmark data. |
| Incomplete validation | Added periodic seam audit test, distinct weighted-sampler tests, and more reference metrics (discrete target, clearance, pair count, mean NND, radius scale). | Execution pending. |
| Extraneous comments / title formatting | Removed two leftover development comment blocks. The publication-style function no longer overwrites the title interpreter. | Source inspected only; visual check in MATLAB pending. |
| No citation/license/release gate | Added `CITATION.cff` and a release checklist. Provided explicit notice that license selection requires owner approval. | The owner must choose a license before public release. |
| Reviewer-specific studies absent | Kept outside the primary benchmark; provenance retained separately. | **Still incomplete** for full-paper reproduction. |

## Scientific cautions

The original power-key selector is **numerically degenerate** for many small probabilities and may cause biased candidate selection. This is a defect in the selection arithmetic, not just a plotting or style concern. The log-key selector implements the intended weighted-without-replacement random ranking but changes which candidates are selected and may change the final geometry and aggregate results. The current reference was produced using the legacy behavior; do not silently switch algorithms in the published reproduction package.

For full-study reproducibility, every experiment associated with the manuscript should have its own validated runnable entry point and its required initial conditions/reference outputs. It would be misleading to describe the present benchmark-only package as a complete reproduction of all reviewer results.

The checkpoint interface uses configuration equality; changing parameters or moving the repository into a new directory will require clearing the old checkpoints or using fresh output directories. A simulation that crashes inside a realization must rerun **that realization**, not resume solver iterations internally.
