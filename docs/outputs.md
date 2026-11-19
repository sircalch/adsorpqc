# Outputs

`adsorpqc assess` and `adsorpqc demo` write these files to the output directory:

| File | Content |
|---|---|
| `report.html` | Self-contained report with every check, its status and diagnostics |
| `adsorpqc_isotherm_fit.svg`, `.pdf`, `.png` | Isotherm data with the fitted models |
| `adsorpqc_gcmc_burnin_trajectory.svg`, `.pdf`, `.png` | Loading per cycle with the burn-in cut-off shaded (when a GCMC series is available) |
| `adsorpqc_summary_table.csv`, `.tex` | Summary table, ready for a spreadsheet or a LaTeX manuscript |
| `methods_snippet.txt` | Draft methods paragraph describing the checks and their results |
| `citation.bib` | BibTeX entry for AdsorpQC |

PNG figures are written at 300 DPI. The methods paragraph is a starting point: edit it so that it describes your study accurately.
