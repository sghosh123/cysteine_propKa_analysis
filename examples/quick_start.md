# Quick-start examples

Run these commands from the repository root after `pip install -e ".[analysis,test]"`.

## 1. Inspect a pKa series

```bash
cpa summarize data/raw/propka_series/Human/Wilsons1/pKa_Cys18P_Wilson_MBD1_apo.dat --ph 7.0
```

The summary distinguishes the observed fraction of framewise *predicted pKa values* below the supplied pH from the separate idealized Henderson–Hasselbalch descriptor.

## 2. Replot a pKa series without legacy hardcoded counts

```bash
cpa plot data/raw/propka_series/Human/Wilsons1/pKa_Cys18P_Wilson_MBD1_apo.dat \
  --output-prefix outputs/wilson_mbd1_cys18
```

This creates `outputs/wilson_mbd1_cys18_trace.png` and `_histogram.png` using all validated input rows.

## 3. Compare all curated traces

```bash
cpa compare --manifest data/manifest.csv --ph 7.0 --output outputs/all_pka_summaries.csv
```

Before interpreting this as a scientific comparison, confirm that residue labels, chain IDs, state definitions, input structures, and PROPKA method/version are compatible across files.

## 4. Parse a static PROPKA report

```bash
cpa parse-report data/raw/propka_reports/Human/Wilsons1/out.pka \
  --output outputs/wilson_mbd1_cysteine_report.csv
```

This extracts primary CYS rows only. It does not infer whether the report represents an apo or copper-bound physical model beyond the input metadata.

## 5. Analyze the local SG environment

```bash
cpa neighbors --topology path/to/system.psf --trajectory path/to/trajectory.dcd \
  --segid P --resid 18 --atom SG --cutoff 4.0 --output outputs/cys18_sg_contacts.csv
```

The output is one row per partner atom per frame within the cutoff. Distances are in Å. A contact within an SG-centered cutoff is not automatically an H-bond or a proton-transfer event.
