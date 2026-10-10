# Source archive audit

## Source archives

- `propKa_cys.zip`: pKa series, PROPKA output reports, VMD Tcl frame-extraction scripts, histogram scripts, and duplicated neighbor-analysis scripts.
- `Wilson_MBD1_apo.zip`, `Wilson_MBD2_apo.zip`, and `Wilson_MBD3_apo.zip`: cysteine-pair distance outputs and MD preparation/production files complementary to the pKa series.

Only framewise pKa data, PROPKA reports, and three distinct cysteine-pair distance traces are copied into this repository. Historical simulation configs and structure preparation files remain outside the maintained runtime path because needed topologies, trajectories, force-field files, and exact run provenance are not all available here.

## Corrections and retained caveats

1. **Hardcoded histogram scripts replaced.** Legacy scripts set a fixed `N` (e.g. 500, 505 or 1005), read exactly that many lines, and embed a dataset-specific file name/title/output path. This can truncate or mislabel data. The new reader uses the actual row count and validates every row.
2. **A stale histogram script was found in the Wilson MBD archives.** `openmm_production_results/histogram.py` refers to a yeast complex Cys13–Cys16 trace rather than the MBD-specific data. It is not migrated; one common plotting implementation is used instead.
3. **A second historical histogram script was broken.** `histogram-version1.py` in the source review included an undefined `ax`/`start`/`end` branch in some variants. Those versions are not copied as executable tools.
4. **Neighbor implementations were duplicated.** Several `find_neighbors.py` copies rebuilt selections in each frame, used a hardcoded empty input path in `main()`, and serialized per-frame atom lists with little provenance. The replacement accepts explicit inputs/selections, resolves the target uniquely, uses the current frame's valid cell when available, includes distances and `resindex`, and exports CSV.
5. **Frame labels and physical time are not interchangeable.** The legacy pKa extraction TCL writes `im + 1`; distance traces in the Wilson archives begin with index 0. No physical time conversion is made without a verified frame stride/timestep mapping.
6. **Dataset length mismatch.** Wilson MBD3 pKa traces have 500 rows, whereas its cysteine-pair distance trace has 505 rows. Do not correlate those series by row number without verifying exact frame correspondence.
7. **Byte-identical MBD2 pair-distance duplicate omitted.** The files `time_vs_dist_500ns_gpu_apo_Wilson_MBD2_cys18_cys21.dat` and its `_version2` copy had identical contents in the source archive. Only one copy is included.
8. **Four short pKa files retained but not assigned a trajectory interpretation.** Four root-level `.dat` files have five rows. They are catalogued as `short_series_unclassified` rather than silently assumed to be MD frame series.
9. **PROPKA version/model caution.** The archived report header states PROPKA 3.0 and warns that it is not identical to PROPKA 2.0 and does not work with ligands. The corresponding VMD scripts export `protein`, omitting Cu and other non-protein atoms. These historical estimates cannot be treated as explicit metal-aware pKa values. Current PROPKA documentation discusses protein–ligand calculations; a newer version still needs model-specific verification for transition-metal sites.
10. **Terminal-value consistency audit.** A script compared the last pKa row in every data series with the matching CYS/chain record in a mirrored `out.pka` report. 40/44 matched to numerical tolerance. Four plain yeast complex series did not: `pKa_Cys13C_yeast_MBD_complex.dat` (+0.73), `pKa_Cys15M_yeast_MBD_complex.dat` (+0.11), `pKa_Cys16C_yeast_MBD_complex.dat` (−1.14), and `pKa_Cys18M_yeast_MBD_complex.dat` (+1.05), where values are series terminal pKa minus report pKa. The `_2` and `_final` variants are consistent with the available report for these residues. This does not identify which version is scientifically correct; those four remain flagged rather than overwritten.
11. **Chain suffix is not a protonation label.** Names such as `Cys18P` and `Cys15M` map to lines such as `CYS 18 P` in PROPKA output, so `P`/`M` is the reported chain/segment token in these files, not a protonation-state notation.
12. **No raw trajectory reproducibility.** The source archives do not contain all referenced topology/coordinate/trajectory files, and no attempt is made to claim a rerun of the original data.

## Data inventory

- 44 original pKa `.dat` files: 40 have at least 200 rows; 4 have five rows.
- 20 original `.pka` report files.
- 3 distinct cysteine-pair distance traces (Wilson MBD1, MBD2, MBD3).

SHA-256 values and original ZIP member names are recorded in the CSV manifests under `data/`. The terminal consistency audit is stored in `data/terminal_frame_pka_report_check.csv` and can be regenerated with `cpa audit-terminal --output outputs/terminal_report_consistency.csv`.
