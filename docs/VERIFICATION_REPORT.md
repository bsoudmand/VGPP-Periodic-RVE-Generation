> Historical packaging/preliminary-check snapshot. Current execution results and release blockers are in `VALIDATION_REPORT.md`.

# VGPP verification report

Date: 30 October 2026. Software version: `0.3.0-preview`.

Repository: https://github.com/bsoudmand/VGPP-Periodic-RVE-Generation

## Repository and publication controls

- The repository was confirmed private in the authenticated GitHub interface.
- No release was published; the repository showed zero tags.
- No root `LICENSE` was added. MIT remains a proposal pending the owner's approval.
- No Zenodo deposit or software DOI was created in this work. No DOI is asserted in the citation metadata.
- The renamed repository retained the original uploaded files and commit history. The preceding documentation commit `d5cf2512ae7d1fd5764fe607d8b3839a1f60b81f` changed only `README.md` and `CITATION.cff`, as verified in the GitHub diff.
- The initial 62-file upload was previously downloaded and every file verified against the source folder by SHA-256. A fresh archive download after renaming was blocked by Chrome (`ERR_BLOCKED_BY_CLIENT` for the GitHub archive host); no security setting was changed or bypassed. This report does not claim a fresh post-rename archive hash verification.

## MATLAB execution

MATLAB executable: MATLAB R2024b executable.

The first batch startup failed before project execution with `failed to load settings errors_warnings plugin`. A second launch used an isolated `MATLAB_PREFDIR` under this task's workspace and completed successfully with exit code 0. Existing user preferences were not modified.

Reported environment:

- MATLAB `24.2.0.2712019 (R2024b)`.
- Java runtime available.
- Image Processing Toolbox license test returned true; `bwdist` and `imregionalmax` were found.
- Active-source fingerprint: `9baf174d6d47a97689dfe9ed97fa7dc655798653aab6d9d1b43a5168bd4016a9`, matching the original manifest.
- `vgpp_self_test(params)` passed.

## Limited numerical smoke test

The unmodified source was executed from a workspace copy. Default model parameters were used with nominal particle area fraction 0.40, target index 1, realization 1, and seed 7401. `run_benchmark_case` was called directly for VGPP and RSA-PBC. This was not execution of the full `main.m` workflow or the README's two-realization/six-target smoke-test configuration.

| Check | VGPP | RSA-PBC |
|---|---:|---:|
| Particle count | 393 | 393 |
| Final radius | 0.18 | 0.18 |
| Achieved area fraction | 0.400025275766895 | 0.400025275766895 |
| Target attained | Yes | Yes |
| Geometry valid | Yes | Yes |
| Reported maximum overlap | 0 | 0 |
| Reported overlap-violating pairs | 0 | 0 |
| Independent minimum clearance | 0.0014457831325292991 | 0.00002833241297323097 |

Independent Python checks evaluated every pair of exported particles with the periodic minimum-image distance in the 10-by-10 domain. Both coordinate sets had 393 particles, centers in the canonical periodic domain, the prescribed radii, and no overlap exceeding `1e-8` model-length units.

The two smoke-test records were matched to the archived 240-record CSV by method, nominal area fraction, realization, and seed. Four discrete fields (`TargetN`, `ReachedTarget`, `NumberOfParticles`, and `OverlappingPairs`) were checked exactly. Eleven continuous fields were checked at absolute tolerance `1e-7`: `FinalVF`, `DiscreteTargetVF`, `MeanRadius`, `RadiusScale`, `MaxOverlap`, `MinClearance`, `LocalVFRSD`, `LowVFCellFraction`, `MeanNND`, `NNDCV`, and `RDFFirstPeak`. Both records passed. The maximum absolute difference was 0 for VGPP and approximately `9.99e-16` for RSA-PBC. Runtime was excluded. This comparison was performed independently in Python; the full `verify_reference_results` routine was not executed against a completed rerun dataset.

## Documentation and metadata review

The README's defaults, seed policy, integer target-count convention, numerical tolerances, toolbox use, output names, and limitations were checked against the source. Version `0.3.0-preview` and both author names were retained. GitHub successfully rendered the citation metadata as an APA software citation with both authors, the software title, version, and correct repository URL. No affiliation, ORCID, release date, journal reference, publication DOI, or software DOI was invented. The public-deposit manuscript and reviewer texts are explicitly conditional on actual public access; their current-status alternatives describe the private repository.

No additional software capability, statistical representativeness, optimum packing density, universal target attainment, or complete manuscript reproducibility is established by these checks.

## Outstanding checks

- Full six-target, 20-realization, two-method benchmark (240 cases) and complete archived-reference comparison.
- Dense-packing paths, stopping/failure behavior, seed-radius growth/recovery, and legacy/log-key comparisons across the study.
- End-to-end `main.m` execution, checkpoint creation/resume/mismatch handling, and figure-generation verification.
- Optional LS studies and supplementary experiments/inputs not fully integrated into the consolidated workflow.
- Owner approval of license, public visibility, and any permanent release or Zenodo publication.

The successful limited tests support initial execution and low-density numerical consistency only. They do not close the full-validation or public-release gates.

MathWorks reference for the isolated preferences directory: https://kr.mathworks.com/matlabcentral/answers/93696-how-do-i-change-the-matlab-preferences-directory-location
