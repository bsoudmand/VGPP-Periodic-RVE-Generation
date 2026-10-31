# Publication materials — published v1.0.0

Manuscript: “Void-guided generation of high-area-fraction periodic particulate RVEs with quantified spatial statistics.” Revised submission to Scientific Reports. No article DOI or publication date is asserted.

Authors: Behzad Hashemi Soudmand; Rasool Mohsenzadeh. Software version: `1.0.0`. The owner approved MIT licensing, public visibility, the versioned release and Zenodo archiving on 2026-10-31. Repository: https://github.com/bsoudmand/VGPP-Periodic-RVE-Generation Public release: https://github.com/bsoudmand/VGPP-Periodic-RVE-Generation/releases/tag/v1.0.0. Release commit: `4d067c6d71c97ae51a5e9743e73b4963afce6c09`. Version-specific DOI: https://doi.org/10.5281/zenodo.23834138. All-versions DOI: https://doi.org/10.5281/zenodo.23834137.

## Repository description

MATLAB implementation of Void-Guided Particle Packing for two-dimensional periodic particulate RVEs, with a fixed-final-radius VGPP/RSA-PBC benchmark, spatial-statistical analysis, manuscript-study scripts and archived reference data.

## License

The approved MIT License is applied in the root LICENSE file, copyright © 2026 Behzad Hashemi Soudmand and Rasool Mohsenzadeh.

## GitHub release description

This implementation generates and analyzes two-dimensional periodic configurations of monodisperse circular particles using Void-Guided Particle Packing, with fixed-final-radius growth/relaxation and RSA-PBC comparison. The reference design comprises six nominal area fractions from 0.40 to 0.70, 20 realizations per target and method, continuous periodic geometry auditing, historical local-cell statistics, nearest-neighbor/RDF descriptors and compatible checkpoints. All 240 core benchmark cases passed the archived-reference comparison across 15 quantities, continuous periodic geometry, independent descriptor recomputation, restart safeguards, and final MATLAB output/dependency checks.

MATLAB R2024b and Image Processing Toolbox are required for packing. Standalone supplementary scripts and archived numerical data cover high-density reachability, sampling/ablation comparisons, LS, grid/gap/domain sensitivity, ordering, finite-chain Monte Carlo and runtime/scaling. Their verification levels are listed individually; not all generators or MC trajectories were rerun. Target attainment is not guaranteed, overlap is assessed at 1e-8, runtimes are environment-dependent, and complete manuscript reproducibility is not claimed.

## Zenodo deposit metadata

**Title:** Void-Guided Particle Packing (VGPP) for Two-Dimensional Periodic Particulate RVEs: MATLAB Implementation and Numerical Reference Data

**Resource type:** Software. **Version:** `1.0.0`. **License:** MIT. **Creators, in order:** Behzad Hashemi Soudmand (given name Behzad; family name Hashemi Soudmand), Rasool Mohsenzadeh (given name Rasool; family name Mohsenzadeh). Affiliations and ORCIDs have not been supplied and are not inferred.

**Abstract:** This software implements Void-Guided Particle Packing for constructing two-dimensional periodic disk configurations intended for particulate representative volume element studies. The MATLAB workflow combines void-guided candidate selection, batch insertion, seed-radius reduction, gradual growth and relaxation, continuous periodic geometric auditing, spatial statistics and RSA-PBC comparison. Historical boundary-side coordinates are retained for local-cell statistics separately from canonical geometric coordinates. The reference benchmark uses six nominal area fractions, 20 seeds per target and method, and restartable per-case checkpoints. All 240 core benchmark cases passed the archived-reference comparison across 15 quantities, continuous periodic geometry, independent descriptor recomputation, restart safeguards, and final MATLAB output/dependency checks. Original standalone manuscript-study scripts and numerical references are included with explicit verification levels and limitations. MATLAB R2024b and Image Processing Toolbox are required for packing. Unavailable continuation/snapshot inputs and incomplete supplementary generator reruns are documented; complete manuscript reproduction is not asserted.

**Keywords:** Void-Guided Particle Packing; VGPP; representative volume element; periodic boundary conditions; particulate microstructures; disk packing; MATLAB; random sequential adsorption; spatial statistics; reproducibility.

**Related source repository:** https://github.com/bsoudmand/VGPP-Periodic-RVE-Generation. **Release:** https://github.com/bsoudmand/VGPP-Periodic-RVE-Generation/releases/tag/v1.0.0. **Release commit:** `4d067c6d71c97ae51a5e9743e73b4963afce6c09`. **Version-specific software DOI:** https://doi.org/10.5281/zenodo.23834138. **All-versions DOI:** https://doi.org/10.5281/zenodo.23834137. The actual Zenodo record title follows CITATION.cff: “Void-Guided Particle Packing (VGPP) for Two-Dimensional Periodic RVE Generation: MATLAB Implementation.”

## Manuscript Code availability

The MATLAB implementation of Void-Guided Particle Packing and the included numerical reference data are publicly available at https://github.com/bsoudmand/VGPP-Periodic-RVE-Generation. Version v1.0.0 (commit 4d067c6d71c97ae51a5e9743e73b4963afce6c09) is archived in Zenodo at https://doi.org/10.5281/zenodo.23834138. The repository documents the validated core benchmark and the verification status and limitations of the supplementary studies.

## Response to Reviewer 2, Comment 11

We thank the reviewer for requesting access to the implementation. The MATLAB source code and included numerical reference data have been deposited in the public repository https://github.com/bsoudmand/VGPP-Periodic-RVE-Generation and archived as version v1.0.0 in Zenodo (https://doi.org/10.5281/zenodo.23834138). The repository provides execution instructions, parameter definitions, numerical tolerances, validation results, and an explicit inventory of supplementary-study verification. The 240-case core benchmark passed the archived-reference comparison and independent coordinate/statistic audits in MATLAB R2024b. Supplementary verification limits are stated explicitly, and complete manuscript reproducibility is not claimed.

These are finalized replacement texts. Editing a manuscript or response document is a separate file operation; no such document edit is asserted here.
