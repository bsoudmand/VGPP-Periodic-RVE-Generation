> Historical preparation snapshot. The current correction and executed results are documented in `VALIDATION_REPORT.md` and `SOURCE_TRACEABILITY.md`.

# MATLAB R2024b verification protocol for preview v0.3

**Environment:** MATLAB R2024b with Image Processing Toolbox and Java enabled. Run from the repository root. Do not use checkpoint directories created by preview v0.2: v0.3 intentionally rejects them.

1. **Self-test:** Add the repository folders as in `main.m`, construct `params=vgpp_config(repoRoot)`, and set `params.sourceFingerprint=vgpp_source_fingerprint(repoRoot)`. Run `vgpp_self_test(params)`. Inspect all assertions, including the finite-geometry and sampler tests (4,000/50,000 candidates).
2. **Smoke test:** In a temporary copy of the repository, set `params.targetVFSweepRuns=2` in `config/vgpp_config.m`; run `main`. Confirm checkpoints and the raw/summary CSVs appear. Rerun `main` with exactly the same code and configuration; confirm the completed checkpoints resume. Plot outputs from this 2-run smoke test are illustrative and not manuscript-validated.
3. **Checkpoint version check:** In another temporary copy, change a harmless source comment in an active MATLAB function (e.g., `src/core/periodic_center_distance.m`). Rerun with existing checkpoints; expect `VGPP:CheckpointMismatch`. Restore the original source before proceeding. Code edits must not be made to any version used for a comparison with archived results.
4. **Full benchmark:** Start with the unchanged source and a clean output directory. Restore `params.targetVFSweepRuns=20`; run `main`. It attempts 20 VGPP and 20 RSA-PBC realizations per target and can take many hours. Check `GeometryValid` and `CompletionStatus` in the raw CSV.
5. **Reference validation:** Confirm `reference_comparison.csv` exists, inspect every `Pass` and `NonfinitePatternMatch` entry, and resolve mismatches. A failed comparison now prevents generation of final figures by `main`. Compare the runtime distribution separately; elapsed runtimes are not bitwise reproducible across machines.
6. **Sampling-variant study:** In a separate copy with fresh output/checkpoint folders, change `weightedSamplingAlgorithm` to `log_keys` and repeat a controlled benchmark. Treat all resulting data as a new numerical variant. Compare target attainment, overlap results, local statistics, runtime distribution, and seed dependence before considering a manuscript-level algorithm change.
7. **Scope gate:** This repository reproduces the fixed-final-radius VGPP/RSA-PBC benchmark only. High-density, LS, ablation, structural-ordering, clearance, and robustness studies need separately validated runnable modules and archival inputs before claims of full-paper computational reproduction.

No MATLAB validation has been performed in the packaging environment. Do not issue a DOI or mark this a final scientific release until the applicable checks pass and the authors approve a license.
