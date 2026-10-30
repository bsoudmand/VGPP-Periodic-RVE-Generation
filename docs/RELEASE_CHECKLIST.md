# Publication release gate

- [ ] Run MATLAB R2024b with Image Processing Toolbox (not available in packaging environment).
- [ ] Execute smoke test and resume from per-realization checkpoints; modify a source file in a disposable copy and confirm checkpoint reuse is blocked.
- [ ] Validate 20-run benchmark counts, success states, geometry, nearest-neighbor/RDF statistics and nonfinite patterns against reference CSV; investigate every deviation.
- [ ] Compare legacy weighted sampler with opt-in log-key sampler in independent runs. Do not conflate numerical repair with exact reproduction.
- [ ] Integrate, test, and document supplementary experiments essential to the paper (high-density tests, positive clearance, LS, ablation, order, sensitivity, scaling).
- [ ] Choose and approve a software license; add canonical LICENSE file prior to public release.
- [ ] Verify author metadata and citation version; issue a Zenodo DOI only after final release.
- [ ] Tag release with MATLAB version, code commit, software requirements, and validated output manifest.
- [ ] Supply exact figure-generation scripts for every published figure, or document supported scope explicitly.

- [ ] Confirm the 4,000-of-50,000 sampler self-test and nonfinite/out-of-domain geometry tests pass in MATLAB.
- [ ] Regenerate docs/MANIFEST.json after any release changes and verify the archived file count.
