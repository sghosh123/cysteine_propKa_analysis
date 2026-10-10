# External workflows

## Export trajectory frames in VMD

1. Edit `export_frames.tcl`: topology/trajectory paths, file types, atom selection, and frame interval.
2. Run it in VMD's Tcl environment.
3. Inspect several generated PDB files, paying particular attention to chain IDs, residue numbering, and whether the intended atoms are included.
4. Install a compatible PROPKA executable (`propka3` or `propka`) and run:

```bash
cpa trajectory-propka pka_frames --executable propka3 --resid 18 --chain P --output-dir propka_run
```

Omit `--resid` to extract all primary CYS pKa records present in each report. Use `--chain` only together with `--resid`. The output CSV carries both the zero-based trajectory frame index and a one-based frame label; neither is physical time without a verified frame-to-time mapping.

Each frame runs in a separate work directory under `propka_run/reports/`. The exact report and stdout/stderr are retained for inspection. PROPKA errors stop the run by default; `--continue-on-error` records failed frames in `errors.log`.

**Metal caveat:** the provided VMD selection defaults to `protein`, which removes Cu and other non-protein atoms. The archived PROPKA 3.0 reports state that this version does not work with ligands. This pipeline therefore does not produce metal-aware cysteine pKa estimates for Cu-bound structures. Do not compare apo and Cu-bound results as a controlled metal effect unless the model explicitly and validly represents metal coordination.
