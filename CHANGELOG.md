# Changelog

## 0.1.0 — Initial curated repository

- Added a validated pKa-series reader, descriptive summaries, optional Henderson–Hasselbalch descriptors, and trace/histogram plotting.
- Added a PROPKA text parser that filters continuation rows and deduplicates primary CYS records repeated in the report summary.
- Added a configurable per-frame PROPKA runner for already-exported PDB snapshots and a VMD frame-export template.
- Added PBC-aware cysteine-SG neighbor/contact extraction with explicit residue identity, partner metadata, and CSV output.
- Curated 44 pKa series, 20 PROPKA reports, and three distinct cysteine-pair distance traces with source metadata and SHA-256 hashes.
- Added terminal-row provenance audit; retained and documented four unresolved mismatches without altering original values.
- Replaced fixed-count, dataset-specific plotting scripts and duplicated neighbor helper scripts with maintained tools.
- Added 16 automated tests and GitHub Actions test workflow.
