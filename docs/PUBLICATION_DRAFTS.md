# Publication documents — drafts pending validation and approval

Manuscript: “Void-guided generation of high-area-fraction periodic particulate RVEs with quantified spatial statistics.” Submitted to Scientific Reports; no publication DOI or publication date is asserted. Authors: Behzad Hashemi Soudmand and Rasool Mohsenzadeh.

Repository: https://github.com/bsoudmand/VGPP-Periodic-RVE-Generation. Preparation version: `0.3.1-preview`. Repository visibility is private. No release or Zenodo record has been published in this preparation. The full validation has not passed; `v1.0.0` is not prepared or tagged.

## Repository description

MATLAB implementation of Void-Guided Particle Packing for two-dimensional periodic particulate RVEs, with VGPP/RSA-PBC benchmarks, original manuscript study scripts, and archived reference data. Full validation is incomplete; see the validation report.

## License — approval required

MIT is proposed. The complete proposed text and copyright line are provided separately as `LICENSE_MIT_DRAFT.txt` for review. No root LICENSE is added and no MIT grant is applied to this repository. Approval of a license does not authorize public visibility, a GitHub release, or Zenodo publication.

## GitHub release description — draft only

This software provides a MATLAB implementation of Void-Guided Particle Packing for two-dimensional periodic disk configurations intended for particulate representative volume element studies. The core workflow compares monodisperse, fixed-final-radius VGPP with RSA-PBC at six nominal particle area fractions, with 20 realizations per target and method. It includes void-guided insertion, radius growth and relaxation, continuous-coordinate periodic auditing, spatial statistics, compatible per-realization checkpoints, and archived-reference comparison. Original supplementary study scripts and archived numerical datasets are supplied separately from the main workflow, with provenance and execution instructions.

MATLAB R2024b and Image Processing Toolbox are required for the packing workflow. Core and safeguard tests, selected archived-data audits, and archived Figure 4 postprocessing have passed. The 240-case rerun remains incomplete and available records show discrepancies in legacy local cell-bin statistics after canonical endpoint normalization. Some supplementary workflows require original solver states or lack full per-snapshot coordinate evidence. Target attainment is not guaranteed, and dense cases can be expensive. This preparation does not claim complete manuscript reproducibility and is not yet suitable for a v1.0.0 release.

## Zenodo metadata — draft, not deposited

**Title:** Void-Guided Particle Packing (VGPP) for Two-Dimensional Periodic Particulate RVEs: MATLAB Implementation and Manuscript Study References

**Resource type:** Software. **Preparation version:** 0.3.1-preview. The final archive version and date must match the actual approved release.

**Creators:** Behzad Hashemi Soudmand (family name: Hashemi Soudmand; given name: Behzad); Rasool Mohsenzadeh (family name: Mohsenzadeh; given name: Rasool). Affiliations and ORCID identifiers are not supplied and are not inferred.

**Abstract:** This software implements Void-Guided Particle Packing in MATLAB for constructing and analyzing two-dimensional periodic configurations of circular particles intended for particulate RVE studies. The consolidated core addresses a monodisperse, fixed-final-radius VGPP/RSA-PBC benchmark using void-guided candidate selection, batch insertion, seed-radius reduction, gradual radius growth, relaxation, and continuous-coordinate periodic geometric audits. The default design comprises six nominal particle area fractions from 0.40 to 0.70, 20 realizations per target and method, and per-realization checkpoints. Original standalone supplementary study families and archived numerical reference data cover high-density packing, comparator and ablation studies, LS benchmarking, gap/grid/domain sensitivity, structural ordering, finite-chain Monte Carlo, and runtime/scaling investigations. MATLAB R2024b and Image Processing Toolbox are required for packing. Targeted MATLAB tests and selected archived-data audits have passed, but full benchmark validation remains incomplete. Canonical endpoint storage produces documented differences in legacy local-bin statistics, and some continuation and per-snapshot evidence remains outside the clean package. Complete manuscript reproducibility is not claimed.

**Keywords:** Void-Guided Particle Packing; VGPP; representative volume element; periodic boundary conditions; particulate microstructures; disk packing; MATLAB; random sequential adsorption; Lubachevsky–Stillinger; spatial statistics; reproducibility.

**Related source repository:** https://github.com/bsoudmand/VGPP-Periodic-RVE-Generation. Record the exact approved commit and release tag when finalizing the archive. **License:** pending owner approval; MIT proposed. **DOI:** none assigned or verified in this preparation. Do not substitute an article DOI for a software DOI.

After approval and validation, use the official GitHub–Zenodo integration where available to archive the exact approved release. Verify that the DOI resolves to the intended version and its archive matches the approved source. Do not publish a permanent archive of this unvalidated preparation without separate authorization.

## Manuscript Code availability — conditional draft

Use the following only after the repository is actually public and its released scope has been verified:

“The MATLAB implementation of Void-Guided Particle Packing and associated numerical reference data are available at https://github.com/bsoudmand/VGPP-Periodic-RVE-Generation. The repository documents the supported benchmark and supplementary study workflows, software requirements, numerical tolerances, validation results, and limitations. The released scope does not establish complete reproduction of every manuscript result.”

After a Zenodo archive is published and verified, add its actual version-specific DOI and released version/commit. No DOI placeholder is presented as an identifier. If validation remains incomplete at publication, state that limitation explicitly in the final manuscript text rather than using the above as an unqualified reproducibility claim.

## Reviewer 2, Comment 11 — conditional draft

Use only after public access is approved and verified:

“We thank the reviewer for requesting access to the implementation. We have made the MATLAB source code and associated numerical reference data publicly accessible at https://github.com/bsoudmand/VGPP-Periodic-RVE-Generation. The repository provides execution instructions, parameter definitions, numerical tolerances, a validation report, and an inventory of the included benchmark and supplementary investigations. The documentation distinguishes verified results from incomplete or unavailable workflows; we do not claim complete manuscript reproducibility beyond the demonstrated released scope.”

Add the actual Zenodo DOI only after archive publication and resolution have been confirmed. State that the manuscript Code availability section has been updated only after the manuscript edit is performed. The repository is currently private, so this response is not ready for submission as a completed-publication statement.
