# Runtime reproducibility investigation

The archived VGPP realization at target 0.60, run 18, seed 7718 records 26295.3338769 s. The corresponding earlier consolidated execution recorded 60.4780564 s. The checked particle counts, area fraction, geometric measures, NN and RDF descriptors matched the archive; this does not prove identical internal iteration histories or original coordinate trajectories.

| Check | Evidence and limit |
|---|---|
| Environment | New executions use MATLAB 24.2.0.2712019 R2024b, Windows 64-bit and Image Processing Toolbox 24.2. The original script is labeled R2024b, but its original hardware, exact build and workload are not recorded for this case. |
| Parameters and seeds | Same 10×10 box, r=0.18, 250×250 grid, nominal target 0.60, integer count convention, seed 7718, legacy sampler, and retained growth/relaxation/stopping settings. The new migration guard checks all numerical configuration fields. |
| Timing scope | Original `code_refined(8).m` places tic/toc directly around `run_fixed_final_radius_VGPP`; NN/RDF/local postprocessing follows the timer. New `run_benchmark_case` also times the generator, before postprocessing and checkpoint writing. New output capture/audit metadata adds small work inside the generator; no evidence attributes the large discrepancy to this overhead. |
| Reuse | Generation timing is measured before checkpoint writing. Resume returns the saved measured time; it does not replace it with the short cache-load time. Migrated v3 records identify their time as preserved from the original generation. |
| Numerical changes | Historical correction adds passive output capture only. Generator source comparison and exact canonical-table checks against saved cases support unchanged generation behavior for checked cases. The earlier strict-audit/endpoint corrections remain documented. |
| Iteration/profiling evidence | Original per-case growth/relaxation iteration trace, profiler output, CPU time, pause history and power-state record are unavailable. No definitive internal effort comparison can be made. |
| Execution conditions | New timings include concurrent audit work; two early launchers briefly overlapped. New seed 7813 at fraction 0.65 records 33819.492411 s. An attempt to inspect Windows power/suspend event logs was denied. Sleep/suspension is a hypothesis, not an established cause of either outlier. |

Conclusion: the specific cause of the 7718 discrepancy remains unresolved. Retain original manuscript timing evidence unchanged. New values are descriptive execution observations, not isolated performance validation or proof of a speedup. Exact wall-clock equality is not an acceptance criterion for geometric/statistical reproduction. Controlled profiling on a recorded environment would be a separately identified timing study; it must not silently replace the original benchmark.

## Corrected-run repetitions

| VGPP seed | Target | Archived runtime (s) | Prior consolidated runtime (s) | Corrected-run runtime (s) |
|---|---:|---:|---:|---:|
| 7718 | 0.60 | 26295.3338769 | 60.4780564 | 61.3857977 |
| 7813 | 0.65 | 132.5114778 | 33819.492411 | 133.0258385 |

The long timing measurements were not repeated in these corrected executions. Checked numerical quantities agree, and the driver verifies canonical configurations against preserved prior checkpoints. This supports reporting the timing outliers as unresolved execution observations; it does not identify their cause or establish a speedup. Original case-level profiling, CPU-time and power-state evidence remains unavailable. New timings include non-isolated workload and some retained migration timings.
