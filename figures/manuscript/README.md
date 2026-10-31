# Author-approved manuscript Figure 4

`Fig4_author_approved.png` is the PNG supplied and approved by the authors in this session. It is preserved unchanged; its checksum is recorded in `FIG4_PROVENANCE.json`. This supersedes the earlier manuscript image used during the unresolved-format investigation. The statistical/display convention is now confirmed by the authors.

Panel (a) reports nearest-neighbor-distance coefficient of variation and the first RDF peak. Panel (b) reports arithmetic mean elapsed time ± sample standard deviation from the original 20-realization benchmark. The displayed runtime range is 0–2000 s. At nominal area fraction 0.60, mean = 1360.393253165 s and sample SD = 5869.0863613567 s; the upper limit is 7229.4796145217 s. The triangular marker identifies upper-whisker clipping. The full numerical SD is retained in `Fig4_original20_statistics.csv`; it is not replaced by the displayed whisker length. RSA-PBC points include attained targets only.

The numerical source is `inputs/reference/Table_R1_fixedRadius_raw.csv` (240 method/target/realization records, including 20 successful VGPP realizations per target). These are the original records, not the expanded 40-realization timing study or newly measured validation runtimes.

The preserved standalone source `studies/manuscript/R5_C13_Regenerate_Fig4_R2024b.m` implements this convention. Its archived-data execution in MATLAB R2024b was verified and its exported statistics matched the original-data table. The core helper `src/benchmark/plot_runtime_fig4.m` implements the same convention. The approved PNG is byte-for-byte identical to the current project export at the source script's expected output location. The original R5-C13 source checksum agrees with the previously recorded provenance, and the packaged version differs only by the documented I/O-prefix and encoding normalization. The previously executed independent rendering has the same dimensions and numerical/display convention but is not pixel-identical. Thus the available plotting source and numerical convention are established; exact raster replay is not claimed. No replacement image has been generated for the approved asset.

For study-level reproduction, follow `docs/MANUSCRIPT_STUDIES.md` to stage the original CSV with the R5-C13 script in a separate writable workspace. Preserve this approved reference image. Figures produced from new benchmark runtimes under `outputs/` are diagnostic outputs and are not substitutes for the manuscript's archived timing figure.

Suggested caption wording for panel (b): “Arithmetic mean runtime ± sample standard deviation from the original 20-realization benchmark; RSA-PBC results include attained targets only. The runtime axis is limited to 0–2000 s. The triangle at area fraction 0.60 indicates that the upper mean-plus-SD limit exceeds the displayed range.”

The revised manuscript PDF was checked: page 11 defines the SD bars and page 12 explains the triangle. The PDF checksum is recorded in the provenance JSON.

## Optional numerical reproduction in a separate output folder

From the repository root in MATLAB R2024b, the existing helper can read the original CSV directly:

```matlab
addpath(fullfile(pwd,'src','benchmark'));
checkDir = fullfile(pwd,'outputs','fig4-archived-data-check');
plot_runtime_fig4(fullfile(pwd,'inputs','reference', ...
    'Table_R1_fixedRadius_raw.csv'), checkDir, checkDir);
```

This creates verification outputs under ignored `outputs/`, leaving the approved reference asset in this folder intact. To use the standalone R5-C13 script instead, copy that script and `inputs/reference/Table_R1_fixedRadius_raw.csv` into a separate writable workspace, then run the copied script. Do not substitute the new benchmark timing CSV for the original-data figure.
