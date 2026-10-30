# Second independent review: technical corrections

This preview addresses the four outstanding findings from the second external review. **It has not been executed in MATLAB R2024b.**

1. **Geometry audit:** `final_overlap_audit_PBC.m` checks finite center coordinates, positive finite radii, matching array lengths, valid domain and tolerance settings, and canonical periodic coordinates before computing any pair distance. Malformed coordinates raise an error and cannot pass the audit. The self-test includes NaN, infinity, boundary and zero-radius cases.
2. **Checkpoint fingerprint:** `vgpp_source_fingerprint.m` calculates SHA-256 from the active benchmark source tree (main/config/core/analysis/benchmark). The fingerprint and format version 2 are stored in each checkpoint and checked on resume. Old preview checkpoints are intentionally incompatible.
3. **Reference comparison and figure order:** `verify_reference_results.m` compares finite/NaN/infinity masks and the signs of infinities. `main.m` now checks a full-run result against the reference CSV before generating any figures. Plot-only instructions were extracted to `vgpp_generate_benchmark_figures.m`; the numerical benchmark computations were not altered by this separation.
4. **Sampler test:** `vgpp_self_test.m` checks uniqueness for 4,000 selected candidates from a population of 50,000 in both the legacy and log-key modes, and checks broad index coverage for log-key sampling. The default legacy calculation has not changed.

## Outstanding release gates

- Execution and quantitative comparison in MATLAB R2024b remain necessary.
- The software owner must approve and supply a LICENSE.
- This archive does not constitute a complete executable reproduction of all reviewer-specific experiments in the paper.
- Because the benchmark source checksum changed, do not resume `.mat` checkpoints from older packages.
