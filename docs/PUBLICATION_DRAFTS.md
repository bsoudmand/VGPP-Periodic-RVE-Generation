# Publication materials — approved v1.0.0 preparation

Manuscript: “Void-guided generation of high-area-fraction periodic particulate RVEs with quantified spatial statistics.” Revised submission to Scientific Reports. No article DOI or publication date is asserted.

Authors: Behzad Hashemi Soudmand; Rasool Mohsenzadeh. Software version: `1.0.0`. The owner approved MIT licensing, public visibility, the versioned release and Zenodo archiving on 2026-10-31. Repository: https://github.com/bsoudmand/VGPP-Periodic-RVE-Generation (publication in progress). Release and archive identifiers will be added after verification.

## Repository description

MATLAB implementation of Void-Guided Particle Packing for two-dimensional periodic particulate RVEs, with a fixed-final-radius VGPP/RSA-PBC benchmark, spatial-statistical analysis, manuscript-study scripts and archived reference data.

## License

The approved MIT License is applied in the root LICENSE file, copyright © 2026 Behzad Hashemi Soudmand and Rasool Mohsenzadeh.

## GitHub release description — draft

This implementation generates and analyzes two-dimensional periodic configurations of monodisperse circular particles using Void-Guided Particle Packing, with fixed-final-radius growth/relaxation and RSA-PBC comparison. The reference design comprises six nominal area fractions from 0.40 to 0.70, 20 realizations per target and method, continuous periodic geometry auditing, historical local-cell statistics, nearest-neighbor/RDF descriptors and compatible checkpoints. All 240 core benchmark cases passed the archived-reference comparison across 15 quantities, continuous periodic geometry, independent descriptor recomputation, restart safeguards, and final MATLAB output/dependency checks.

MATLAB R2024b and Image Processing Toolbox are required for packing. Standalone supplementary scripts and archived numerical data cover high-density reachability, sampling/ablation comparisons, LS, grid/gap/domain sensitivity, ordering, finite-chain Monte Carlo and runtime/scaling. Their verification levels are listed individually; not all generators or MC trajectories were rerun. Target attainment is not guaranteed, overlap is assessed at 1e-8, runtimes are environment-dependent, and complete manuscript reproducibility is not claimed.

## Zenodo deposit metadata — draft, not deposited

**Title:** Void-Guided Particle Packing (VGPP) for Two-Dimensional Periodic Particulate RVEs: MATLAB Implementation and Numerical Reference Data

**Resource type:** Software. **Version:** `1.0.0`. **License:** MIT. **Creators, in order:** Behzad Hashemi Soudmand (given name Behzad; family name Hashemi Soudmand), Rasool Mohsenzadeh (given name Rasool; family name Mohsenzadeh). Affiliations and ORCIDs have not been supplied and are not inferred.

**Abstract:** This software implements Void-Guided Particle Packing for constructing two-dimensional periodic disk configurations intended for particulate representative volume element studies. The MATLAB workflow combines void-guided candidate selection, batch insertion, seed-radius reduction, gradual growth and relaxation, continuous periodic geometric auditing, spatial statistics and RSA-PBC comparison. Historical boundary-side coordinates are retained for local-cell statistics separately from canonical geometric coordinates. The reference benchmark uses six nominal area fractions, 20 seeds per target and method, and restartable per-case checkpoints. All 240 core benchmark cases passed the archived-reference comparison across 15 quantities, continuous periodic geometry, independent descriptor recomputation, restart safeguards, and final MATLAB output/dependency checks. Original standalone manuscript-study scripts and numerical references are included with explicit verification levels and limitations. MATLAB R2024b and Image Processing Toolbox are required for packing. Unavailable continuation/snapshot inputs and incomplete supplementary generator reruns are documented; complete manuscript reproduction is not asserted.

**Keywords:** Void-Guided Particle Packing; VGPP; representative volume element; periodic boundary conditions; particulate microstructures; disk packing; MATLAB; random sequential adsorption; spatial statistics; reproducibility.

**Related source repository:** https://github.com/bsoudmand/VGPP-Periodic-RVE-Generation. Add the actual approved release URL and exact commit only after verification. **Software DOI:** none created or verified. Do not substitute an article DOI. Follow `ZENODO_ARCHIVING_PLAN.md` after validation and authorization.

## Manuscript Code availability — conditional publication template

Do not submit this template while the repository is private. Replace the bracketed fields only after actual publication and identifier verification:

“The MATLAB implementation of Void-Guided Particle Packing and the included numerical reference data are available at https://github.com/bsoudmand/VGPP-Periodic-RVE-Generation (version [approved release tag], commit [verified release commit]). The versioned software archive is preserved in Zenodo ([verified version-specific DOI]). The repository documents the validated core benchmark and the verification status and limitations of the supplementary studies.”

## Reviewer 2, Comment 11 — conditional publication template

Use only after public access and the archive are verified:

“We thank the reviewer for requesting access to the implementation. The MATLAB source code and included numerical reference data have been deposited in the public repository https://github.com/bsoudmand/VGPP-Periodic-RVE-Generation (version [approved release tag]) and archived in Zenodo ([verified version-specific DOI]). The repository provides execution instructions, parameter definitions, numerical tolerances, validation results, and an explicit inventory of supplementary-study verification. The reproducibility statement is limited to the workflows supported by the reported validation.”

State that the manuscript Code availability section was updated only after that edit is actually performed. No publication DOI, public-release URL, software DOI or validation outcome is invented in these templates.
