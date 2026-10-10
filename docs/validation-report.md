# Validation report for repository snapshot 0.1.0

## Automated checks performed

- **16 tests passed** under the development environment.
- All 44 curated pKa series passed the strict two-column parser and matched the row counts in `data/manifest.csv`.
- 20 archived `.pka` report files were preserved with source-member and SHA-256 metadata.
- The PROPKA text parser was checked against the archived report layout, including repeated detailed/summary records and continuation rows.
- The local-environment module was exercised with a synthetic MDAnalysis Universe containing a boundary-crossing contact; the minimum-image distance was recovered as expected and hydrogens were excluded by default.
- The per-frame PROPKA orchestrator was exercised with a mock executable. No real PROPKA calculations were launched.
- `cpa compare --manifest data/manifest.csv --ph 7.0` successfully summarized all 44 time series.
- `python -m compileall -q src tests` completed successfully.
- A package wheel was built successfully with `python -m pip wheel . --no-deps --no-build-isolation`.
- `cpa --help` successfully listed all seven available subcommands.
- `cpa audit-terminal` checked the terminal pKa values against mirrored reports: 40 matched and 4 were flagged as mismatches.
- `cpa parse-report` parsed the real archived Wilson MBD1 PROPKA report and returned four unique cysteine residue/chain records.
- Example trace and histogram plots were generated with the unified plotting code.

## Unresolved checks

- No real VMD session was available to execute `scripts/export_frames.tcl`.
- No real PROPKA binary was invoked on exported snapshots.
- Full source topologies and MD trajectories were not available to reproduce the historical time-series data.
- The four provenance mismatches remain unresolved pending source-run history; the series were not edited.
- Distance-trace units and frame correspondence are recorded as assumptions where the source workflow does not fully establish them.
- No license has been chosen; the snapshot should not be presented as an openly licensed public project until that decision and result-data redistribution permissions are resolved.
