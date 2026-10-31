# Author-approved Figure 4 provenance

The authors have supplied and approved the corrected Figure 4 with arithmetic mean ± sample SD and the upper-truncation triangle. The previously reported scientific/display-convention issue is resolved. The supplied PNG is preserved unchanged under `figures/manuscript/Fig4_author_approved.png`; it matches the current project export byte-for-byte. No different manuscript figure has been generated or substituted.

The preserved plotting source is `studies/manuscript/R5_C13_Regenerate_Fig4_R2024b.m`. Its original source checksum agrees with the previously inventoried MATLAB script; the packaged version changes only the machine-specific I/O prefix and encoding/line endings. `src/benchmark/plot_runtime_fig4.m` implements the same numerical/display convention. The R4-C10 logarithmic plot is a historical alternative and is not the approved reference.

Inputs are the unchanged `inputs/reference/Table_R1_fixedRadius_raw.csv`, containing 240 original method/target/realization records and 20 successful VGPP realizations per target. The expanded 40-realization analysis and newly measured validation timings are separate datasets. RSA-PBC points include attained targets only.

At nominal area fraction 0.60, arithmetic mean runtime = 1360.393253165 s, sample SD = 5869.0863613567 s, and the upper mean-plus-SD limit = 7229.4796145217 s. The displayed axis is 0–2000 s. The triangle marks upper-whisker clipping; the exported numerical SD remains unmodified. The maximum individual runtime, 26295.3338769 s, is a different quantity from the plotted mean or SD.

The standalone R5-C13 script was executed in MATLAB R2024b on the original-data copy during this preparation. Independent Python recalculation checked 72 exported counts/statistics across six target fractions, with maximum absolute difference 5.12e-13; all pass 1e-7. Evidence is in `FIG4_STATISTICS_CHECK.csv` and `figures/manuscript/Fig4_original20_statistics.csv`. The revised manuscript PDF was inspected: page 11 defines SD bars and page 12 explains the truncation triangle.

Checksums of the approved PNG, original data, original/packaged source, and revised manuscript PDF are recorded in `figures/manuscript/FIG4_PROVENANCE.json`. The previously executed independent raster is not pixel-identical to the supplied PNG. This does not reopen the author-approved formatting choice; it limits any claim of exact raster replay. The available source and its numerical convention are preserved without changing the approved image.

See `figures/manuscript/README.md` for reproduction/data-staging instructions. New-run plots under ignored `outputs/` are diagnostic benchmark outputs and do not replace the manuscript reference figure.
