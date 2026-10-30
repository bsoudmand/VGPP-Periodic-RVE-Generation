# Audit decisions and unresolved items

- Inventory: 8 separately uploaded MATLAB files and 34 MATLAB files from `R2_C38_RuntimeScaling.zip`. Hashes and dispositions are in `SCRIPT_AUDIT.csv`.
- Many 5,500–5,700-line reviewer studies repeat the same core implementation byte-for-byte or with minor textual differences. The published **original 20-realization benchmark** in `code_refined(8).m` is the basis of the modular entry point. This is a provenance choice, not proof that it reproduces all later experiments.
- Reviewer-specific modified kernels are intentionally **not silently unified**. Some change the particle-size model, convergence settings, relaxation or output protocol. Combining them without rerunning their validated test cases could change results.
- `code_comparison(1).m` implements an earlier CL-ASAA/PBC study with `minGap=0.02`, a 600×600 grid and different stopping rules. It is not the fixed-final-radius VGPP manuscript benchmark.
- `code_Fig_Replot(1).m` uses *volume-fraction* labels and the original figure-number sequence. It is archival. The original manuscript's figure 4b cannot be reconstructed reliably from a visual screenshot alone; plotting must use CSV statistics and disclose clipping conventions.
- Both `code_GA(1).m` and `VGPP_Graphical_Abstract_Final_NoWarnings(2).m` are duplicates. The graphical abstract is an illustration, not proof of geometric validity.
- The supplied 240-row `Table_R1_fixedRadius_raw(1).csv` is packaged as an immutable reference. It supports post-run comparison, but no claim of rerun equivalence is made.
- The supplied ZIP includes no complete set of coordinate snapshots, LS checkpoints, or all reviewer-study input files; these are prerequisites for certain independent audits and exact reproductions.
- Real-world run prerequisite: MATLAB R2024b with Image Processing Toolbox. The user should run the 2-realization smoke test, then the full 20-realization benchmark, and compare against reference CSV. This has not been performed in this environment.
- Known performance consideration: the PBC pairwise geometric operations and collective relaxation at dense packing are computationally expensive. Optimization may change the order of floating-point operations; preserve them unless independently validated.
- The zero/nonzero overlap rule uses a permitted numerical overlap of `1e-8`. This is *not* a guarantee of strictly nonintersecting disks under exact arithmetic.

## Duplicate-kernel audit

`KERNEL_VARIANTS.csv` compares 5 key numerical functions across every supplied script after whitespace normalization. Implementations classified as `Modified` are **not** substituted into the core because equivalence has not been demonstrated.

From the supplied original 20-realization CSV (no code rerun), at nominal phi_A=0.60: mean=1360.393253165 s, sample SD=5869.086361357 s, median=47.604056250 s, maximum individual time=26295.333876900 s. These quantities must not be confused in Fig. 4b.
