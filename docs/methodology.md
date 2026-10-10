# Methodology and interpretation

## 1. Framewise predicted pKa values

The legacy `.dat` series have two whitespace-delimited columns: a frame label and a predicted pKa. The reader rejects malformed/non-finite values, repeated labels, and labels that do not increase strictly. It does not silently infer the physical time step.

Descriptive summaries report the number of rows, first/last labels, mean, median, sample standard deviation, quantiles, extrema, and lag-1 autocorrelation. These quantities characterize the supplied series only. They do not establish equilibrium or convergence, and frame correlation means the number of rows is not the number of independent samples.

### Henderson–Hasselbalch descriptor

For the thiol/thiolate equilibrium \(\mathrm{RSH \rightleftharpoons RS^- + H^+}\), the idealized deprotonated fraction is

\[
  f_{\mathrm{RS^-}} = \frac{1}{1 + 10^{pK_a-pH}}.
\]

The toolkit applies this mapping to each input pKa value only when the user explicitly supplies a pH. The resulting average is named `mean_HH_thiolate_fraction_model_only`. It is not asserted to be the actual protonation fraction in the simulated ensemble: a predicted pKa has model error; the snapshots are correlated; and the pKa method may omit relevant coordination chemistry.

## 2. PROPKA report parsing and reruns

The report parser extracts a primary line of the form `CYS <resid> <chain> <pKa> ...`, excluding continuation lines without a pKa and deduplicating repeated entries in the detailed and final summary tables by residue number and chain. In this archive's file names, suffixes such as `Cys18P` correspond to residue 18 and chain/segment token `P`; the suffix is not a protonation-state label. For new trajectories, `scripts/export_frames.tcl` writes PDB snapshots and a frame manifest. The `cpa trajectory-propka` command runs an installed executable separately in each frame's own work directory, saves stdout/stderr and the `.pka` output, and extracts CYS records. The residue number and PDB chain ID must be checked against the actual exported PDB.

The code cannot ensure that an external pKa method models a copper site correctly. A report in the source archive is from PROPKA 3.0 and includes the warning “does not work with ligands.” Current PROPKA documentation describes protein–ligand prediction, but a ligand capability statement is not proof that a transition-metal coordination environment is parameterized for the user's system. For copper-bound snapshots, confirm the current PROPKA version's handling of the exact Cu model and inspect how the PDB was prepared. If the exported selection excludes Cu, the estimate is not metal-aware.

## 3. Local cysteine environment

`cpa neighbors` uses MDAnalysis `capped_distance` to enumerate partner atoms within a user-set cutoff around one uniquely selected atom (default `SG`). If valid unit-cell dimensions exist, the current frame's `[lx, ly, lz, alpha, beta, gamma]` box is passed so the minimum-image convention can be applied. Distances are reported in Å under MDAnalysis' usual coordinate convention. The code falls back to non-periodic distance calculations when box dimensions are missing or invalid; this is not hidden, so users should inspect the source trajectory's cell metadata.

Hydrogens are excluded by default using element information with an atom-name fallback. For H-bond analysis, an SG-neighbor cutoff alone must not be described as a hydrogen bond. Donor/acceptor identity, hydrogen geometry, and an angular criterion would be needed for a geometry-based H-bond assignment.

## 4. Cysteine-pair distance companion data

The Wilson MBD distance series are included because cysteine–cysteine geometry may provide useful structural context. Their frame index starts at zero, while the PROPKA Tcl output convention is generally one-based. The Wilson MBD3 pair-distance trace has 505 rows while the corresponding curated pKa traces have 500 rows; therefore it must not be joined by row position without checking saved-frame correspondence. The data manifest records this limitation. Pair distance units are presumed to be Å from the original distance workflow, but users should verify the original script and trajectory reader before quantitatively interpreting the values.

## 5. Provenance consistency check

`cpa audit-terminal` compares the terminal pKa in each `.dat` series with the same CYS/residue/chain entry in the `out.pka` report under the mirrored folder. In the included snapshot, 40/44 series match to numerical tolerance. Four plain yeast complex series do not match the report, while corresponding `_2`/`_final` variants match. This check only identifies file-to-file inconsistency: it does not establish which rerun is correct. The unmatched series are retained verbatim and identified in `data/terminal_frame_pka_report_check.csv` rather than silently overwritten.

## 6. Scientific claims

Use descriptions such as “framewise PROPKA-estimated cysteine pKa distribution” and “SG-centered contact occupancy.” Avoid claiming direct proton transfer, an equilibrium protonation population, free-energy differences, or converged pKa values from these summaries alone. Follow-up methods may include independent structure preparation, comparison against a suitable experimental/reference estimate, sensitivity analysis, block-based uncertainty assessment, and a pKa approach that explicitly treats the metal environment.

## References and tool documentation

- PROPKA command documentation: https://propka.readthedocs.io/en/latest/command.html
- PROPKA quick start and output behavior: https://propka.readthedocs.io/en/latest/quickstart.html
- MDAnalysis distance routines and minimum-image convention: https://docs.mdanalysis.org/stable/documentation_pages/lib/distances.html
