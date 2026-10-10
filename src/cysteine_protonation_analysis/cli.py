"""Command-line interface for cysteine-protonation-analysis."""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import sys

from .io import read_pka_series, parse_propka_cysteines
from .metrics import summarize_series
from .plotting import plot_series


def _write_dicts(rows: list[dict[str, object]], path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError("No rows to write")
    fields = list(dict.fromkeys(key for row in rows for key in row.keys()))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _cmd_summarize(args) -> int:
    series = read_pka_series(args.input)
    row = summarize_series(series, ph=args.ph)
    if args.output:
        _write_dicts([row], args.output)
        print(f"Wrote {args.output}")
    else:
        print(json.dumps(row, indent=2))
    return 0


def _cmd_plot(args) -> int:
    series = read_pka_series(args.input)
    trace, hist = plot_series(series, args.output_prefix, bins=args.bins)
    print(f"Wrote {trace}\nWrote {hist}")
    return 0


def _expand_inputs(paths: list[str], manifest: str | None) -> list[Path]:
    if manifest:
        manifest_path = Path(manifest).resolve()
        with manifest_path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        # Repository manifests use project-root-relative paths. For a custom
        # manifest outside a project, paths are relative to the manifest folder.
        base = manifest_path.parent
        for parent in [manifest_path.parent, *manifest_path.parents]:
            if (parent / "pyproject.toml").exists():
                base = parent
                break
        output = []
        for row in rows:
            p = Path(row["path"])
            output.append(p if p.is_absolute() else (base / p).resolve())
        return output
    return [Path(p) for p in paths]


def _cmd_compare(args) -> int:
    inputs = _expand_inputs(args.inputs, args.manifest)
    if not inputs:
        raise ValueError("Supply one or more .dat paths or --manifest")
    rows = [summarize_series(read_pka_series(path), ph=args.ph) for path in inputs]
    _write_dicts(rows, args.output)
    print(f"Summarized {len(rows)} series to {args.output}")
    return 0


def _cmd_parse_report(args) -> int:
    rows = parse_propka_cysteines(args.input)
    if not rows:
        raise ValueError(f"No primary CYS pKa entries found in {args.input}")
    _write_dicts(rows, args.output)
    print(f"Extracted {len(rows)} CYS pKa entries to {args.output}")
    return 0


def _cmd_neighbors(args) -> int:
    try:
        import MDAnalysis as mda
    except ImportError as exc:
        raise RuntimeError("Trajectory neighbor analysis requires MDAnalysis; install with `pip install .[analysis]`") from exc
    from .neighbors import cysteine_atom_contacts, write_contacts_csv
    universe = mda.Universe(args.topology, args.trajectory)
    rows = cysteine_atom_contacts(universe, resid=args.resid, atom_name=args.atom,
                                  segid=args.segid, chain=args.chain, cutoff=args.cutoff,
                                  include_hydrogens=args.include_hydrogens)
    write_contacts_csv(rows, args.output)
    print(f"Wrote {len(rows)} atom contacts to {args.output}")
    return 0


def _cmd_audit_terminal(args) -> int:
    from .provenance import audit_terminal_values
    rows = audit_terminal_values(args.series_root, args.reports_root, tolerance=args.tolerance)
    _write_dicts(rows, args.output)
    counts = {}
    for row in rows:
        counts[str(row.get("status", "unknown"))] = counts.get(str(row.get("status", "unknown")), 0) + 1
    print(f"Audited {len(rows)} series: {counts}; wrote {args.output}")
    return 0


def _cmd_trajectory_propka(args) -> int:
    from .propka_runner import run_propka_frames
    output = run_propka_frames(args.frame_dir, args.output_dir, executable=args.executable,
                               resid=args.resid, chain=args.chain,
                               continue_on_error=args.continue_on_error)
    print(f"Wrote pKa time series: {output}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="cpa", description="Analyze cysteine pKa estimates and local environments.")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("summarize", help="summarize one frame-labelled pKa series")
    p.add_argument("input", help="two-column whitespace-delimited file: frame_label pKa")
    p.add_argument("--ph", type=float, default=None, help="optional pH for a model-only Henderson-Hasselbalch mapping")
    p.add_argument("--output", help="write summary as CSV rather than JSON to stdout")
    p.set_defaults(func=_cmd_summarize)

    p = sub.add_parser("plot", help="plot pKa trace and histogram")
    p.add_argument("input")
    p.add_argument("--output-prefix", default="pka", help="prefix for _trace.png and _histogram.png")
    p.add_argument("--bins", type=int, default=25)
    p.set_defaults(func=_cmd_plot)

    p = sub.add_parser("compare", help="summarize multiple pKa series to one CSV")
    p.add_argument("inputs", nargs="*", help="pKa series paths; mutually exclusive with --manifest")
    p.add_argument("--manifest", help="CSV dataset manifest with a path column")
    p.add_argument("--ph", type=float, default=None)
    p.add_argument("--output", required=True)
    p.set_defaults(func=_cmd_compare)

    p = sub.add_parser("parse-report", help="extract primary CYS pKa lines from a PROPKA .pka text report")
    p.add_argument("input")
    p.add_argument("--output", required=True)
    p.set_defaults(func=_cmd_parse_report)

    p = sub.add_parser("audit-terminal", help="compare series terminal values with mirrored out.pka reports")
    p.add_argument("--series-root", default="data/raw/propka_series")
    p.add_argument("--reports-root", default="data/raw/propka_reports")
    p.add_argument("--tolerance", type=float, default=1e-6)
    p.add_argument("--output", required=True)
    p.set_defaults(func=_cmd_audit_terminal)

    p = sub.add_parser("trajectory-propka", help="run PROPKA on exported frame_XXXXXX.pdb files")
    p.add_argument("frame_dir", help="directory of exported PDB frames")
    p.add_argument("--executable", default="propka3", help="PROPKA executable available on PATH")
    p.add_argument("--resid", type=int, default=None, help="only extract this cysteine residue number")
    p.add_argument("--chain", default=None, help="chain ID to disambiguate --resid")
    p.add_argument("--output-dir", required=True, help="new or existing output directory")
    p.add_argument("--continue-on-error", action="store_true", help="record failed frames and continue")
    p.set_defaults(func=_cmd_trajectory_propka)

    p = sub.add_parser("neighbors", help="extract PBC-aware atom contacts around one cysteine atom")
    p.add_argument("--topology", required=True)
    p.add_argument("--trajectory", required=True)
    p.add_argument("--segid", default=None, help="segment ID; useful for PSF/NAMD-style topology files")
    p.add_argument("--chain", default=None, help="PDB chain ID if available in the topology")
    p.add_argument("--resid", required=True, type=int)
    p.add_argument("--atom", default="SG")
    p.add_argument("--cutoff", type=float, default=4.0, help="cutoff in Å")
    p.add_argument("--include-hydrogens", action="store_true")
    p.add_argument("--output", required=True)
    p.set_defaults(func=_cmd_neighbors)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "compare" and args.manifest and args.inputs:
        parser.error("Use either positional inputs or --manifest, not both")
    try:
        return args.func(args)
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
