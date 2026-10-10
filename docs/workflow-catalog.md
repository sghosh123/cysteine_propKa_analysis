# Workflow catalog

| Workflow | Command | Key assumption | Status |
|---|---|---|---|
| Single pKa summary | `cpa summarize series.dat --ph 7.0` | First column is an ordered frame label; not time | Tested |
| Multi-series summary | `cpa compare ... --output summary.csv` | Compared files share the same meaning of pKa and metadata is not silently harmonized | Tested |
| pKa plots | `cpa plot series.dat --output-prefix outputs/name` | All actual rows are loaded; bin count is configurable | Tested |
| PROPKA text extraction | `cpa parse-report report.pka --output cys.csv` | Extracts primary CYS lines, excluding continuation lines and duplicate summary-table entries | Tested on archived report format |
| Terminal provenance audit | `cpa audit-terminal --output audit.csv` | Assumes mirrored report/data directories and infers target chain/residue from the file name | Tested against all 44 archived series; 4 mismatches flagged |
| Frame snapshot export | `scripts/export_frames.tcl` | User selects valid topology, trajectory, and atom selection in VMD | Syntax reviewed; not run in VMD here |
| Per-frame PROPKA runner | `cpa trajectory-propka ...` | Installed executable and exported PDBs are compatible; chain/residue mapping is correct | Tested with mock executable; real PROPKA not run here |
| Cysteine SG contacts | `cpa neighbors ...` | MDAnalysis-compatible topology/trajectory; valid box for PBC; target atom resolves uniquely | Tested on synthetic Universe, including boundary crossing |
| Cysteine-pair geometry context | `data/reference/cysteine_pair_distances/` | Distances presumed Å; frame mapping not validated for direct correlation | Data included; provenance caveat |

Potential future analysis that should not be performed automatically: correlate pKa and cysteine-pair distances. It requires verified frame mapping, trajectory stride, consistent sampled frames, and a plan for serial correlation before statistical interpretation.
