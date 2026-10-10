"""Cross-check time-series terminal values against colocated PROPKA reports."""
from __future__ import annotations
from pathlib import Path
import re
from .io import read_pka_series, parse_propka_cysteines

_TARGET_RE = re.compile(r"Cys(\d+)([A-Za-z])_", re.IGNORECASE)


def audit_terminal_values(series_root: str | Path, reports_root: str | Path,
                          tolerance: float = 1e-6) -> list[dict[str, object]]:
    """Compare each series' final value to its sibling ``out.pka`` primary row.

    The directory trees under ``series_root`` and ``reports_root`` are expected
    to mirror each other. This only establishes numerical consistency between
    supplied files; it does not prove the report is the source of the series.
    """
    series_root = Path(series_root).resolve()
    reports_root = Path(reports_root).resolve()
    if tolerance < 0:
        raise ValueError("tolerance must be non-negative")
    if not series_root.is_dir() or not reports_root.is_dir():
        raise FileNotFoundError("series_root and reports_root must both be directories")
    results: list[dict[str, object]] = []
    for series_path in sorted(series_root.rglob("*.dat")):
        match = _TARGET_RE.search(series_path.name)
        try:
            series = read_pka_series(series_path)
        except ValueError as exc:
            results.append({"series": str(series_path.relative_to(series_root)), "status": "invalid_series", "detail": str(exc)})
            continue
        relative_parent = series_path.parent.relative_to(series_root)
        report_path = reports_root / relative_parent / "out.pka"
        row: dict[str, object] = {
            "series": str(series_path.relative_to(series_root)),
            "report": str(report_path.relative_to(reports_root)),
            "n_rows": int(series.values.size),
            "last_frame_label": int(series.frames[-1]),
            "last_pka": float(series.values[-1]),
            "target_resid_from_filename": "",
            "target_chain_from_filename": "",
            "report_pka": "",
            "difference_last_minus_report": "",
            "status": "",
            "detail": "",
        }
        if not match:
            row.update(status="target_not_parsed", detail="Could not infer Cys residue/chain from filename")
        elif not report_path.exists():
            row.update(status="report_missing", detail="No matching out.pka report in mirrored directory")
        else:
            resid, chain = int(match.group(1)), match.group(2)
            row["target_resid_from_filename"] = resid
            row["target_chain_from_filename"] = chain
            report_rows = parse_propka_cysteines(report_path)
            matched = [r for r in report_rows if int(r["resid"]) == resid and str(r["chain"]).upper() == chain.upper()]
            if len(matched) != 1:
                row.update(status="target_missing_or_ambiguous", detail=f"Found {len(matched)} matching CYS record(s) in report")
            else:
                report_pka = float(matched[0]["pka"])
                difference = float(series.values[-1] - report_pka)
                row["report_pka"] = report_pka
                row["difference_last_minus_report"] = difference
                row["status"] = "match" if abs(difference) <= tolerance else "mismatch"
                if row["status"] == "mismatch":
                    row["detail"] = "Terminal series value differs from colocated out.pka; do not assume same run provenance"
        results.append(row)
    return results
