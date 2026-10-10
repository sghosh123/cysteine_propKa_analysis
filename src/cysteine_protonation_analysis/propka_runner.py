"""Run an installed PROPKA executable on exported trajectory-frame PDBs.

This runner intentionally does not build or protonate metal-coordination models.
It runs the external executable on each PDB as written and records the exact
report used for extraction. Inspect the structure and PROPKA limitations first.
"""
from __future__ import annotations
from pathlib import Path
import csv
import hashlib
import json
import shutil
import subprocess
import re
from .io import parse_propka_cysteines

_FRAME_RE = re.compile(r"frame[_-]?(\d+)", re.IGNORECASE)


def _write_rows(rows: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["frame_index_zero_based", "frame_label_one_based", "resname", "resid", "chain", "pka", "pka_report"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def run_propka_frames(frame_dir: str | Path, output_dir: str | Path, executable: str = "propka3",
                      resid: int | None = None, chain: str | None = None,
                      continue_on_error: bool = False) -> Path:
    """Run PROPKA per PDB in ``frame_dir`` and write selected CYS predictions.

    Frame PDB names must contain an integer like ``frame_000123.pdb``; that
    integer is interpreted as a zero-based trajectory frame index. Each frame is
    run in an isolated directory, preventing output-file collisions. If ``resid``
    is supplied, only that CYS is returned; a ``chain`` may disambiguate repeats.
    """
    frame_dir = Path(frame_dir).resolve()
    output_dir = Path(output_dir).resolve()
    if not frame_dir.is_dir():
        raise FileNotFoundError(f"Frame directory does not exist: {frame_dir}")
    pdb_files = sorted(frame_dir.glob("*.pdb"))
    if not pdb_files:
        raise ValueError(f"No .pdb frames found in {frame_dir}")
    if resid is None and chain is not None:
        raise ValueError("--chain requires --resid")
    output_dir.mkdir(parents=True, exist_ok=True)
    report_root = output_dir / "reports"
    report_root.mkdir(exist_ok=True)
    rows: list[dict[str, object]] = []
    errors: list[str] = []

    for pdb_path in pdb_files:
        match = _FRAME_RE.search(pdb_path.stem)
        if not match:
            message = f"Cannot infer frame index from filename {pdb_path.name}; use frame_000000.pdb style names"
            if continue_on_error:
                errors.append(message); continue
            raise ValueError(message)
        frame_index = int(match.group(1))
        work_dir = report_root / pdb_path.stem
        work_dir.mkdir(parents=True, exist_ok=True)
        local_pdb = work_dir / pdb_path.name
        shutil.copy2(pdb_path, local_pdb)
        # Clear stale outputs before running so a failed or changed rerun cannot
        # accidentally parse an old pKa report for the current snapshot.
        for old_output in [*work_dir.glob("*.pka"), work_dir / "stdout.txt", work_dir / "stderr.txt"]:
            if old_output.exists():
                old_output.unlink()
        try:
            proc = subprocess.run([executable, local_pdb.name], cwd=work_dir, text=True,
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
            (work_dir / "stdout.txt").write_text(proc.stdout, encoding="utf-8")
            (work_dir / "stderr.txt").write_text(proc.stderr, encoding="utf-8")
            if proc.returncode != 0:
                raise RuntimeError(f"{executable} returned exit status {proc.returncode}: {proc.stderr[-500:]}")
            candidate_reports = sorted(p for p in work_dir.glob("*.pka") if p.name != "stdout_report.pka")
            exact_report = work_dir / f"{local_pdb.stem}.pka"
            report_path = exact_report if exact_report.exists() else (candidate_reports[0] if candidate_reports else None)
            if report_path is not None:
                cys_rows = parse_propka_cysteines(report_path)
                report_label = str(report_path.relative_to(output_dir))
            else:
                stdout_report = work_dir / "stdout_report.pka"
                stdout_report.write_text(proc.stdout, encoding="utf-8")
                cys_rows = parse_propka_cysteines(stdout_report)
                report_label = str(stdout_report.relative_to(output_dir))
            metadata = {
                "source_pdb": str(pdb_path),
                "source_pdb_sha256": hashlib.sha256(pdb_path.read_bytes()).hexdigest(),
                "command": [executable, local_pdb.name],
                "return_code": int(proc.returncode),
                "stdout_file": str((work_dir / "stdout.txt").relative_to(output_dir)),
                "stderr_file": str((work_dir / "stderr.txt").relative_to(output_dir)),
                "report_file": str(report_path.relative_to(output_dir)) if report_path is not None else report_label,
            }
            (work_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
            if resid is not None:
                cys_rows = [r for r in cys_rows if int(r["resid"]) == int(resid)]
            if chain is not None:
                cys_rows = [r for r in cys_rows if str(r["chain"]) == str(chain)]
            if not cys_rows:
                raise ValueError(f"No matching primary CYS pKa line in report for {pdb_path.name}; check residue numbering and PDB chain IDs")
            if resid is not None and chain is None and len({str(r["chain"]) for r in cys_rows}) > 1:
                raise ValueError(f"Residue {resid} appears in multiple chains in {pdb_path.name}; supply --chain")
            for item in cys_rows:
                rows.append({
                    "frame_index_zero_based": frame_index,
                    "frame_label_one_based": frame_index + 1,
                    "resname": item["resname"],
                    "resid": item["resid"],
                    "chain": item["chain"],
                    "pka": item["pka"],
                    "pka_report": report_label,
                })
        except Exception as exc:
            if not continue_on_error:
                raise
            errors.append(f"{pdb_path.name}: {exc}")
    out_csv = output_dir / "cysteine_pka_by_frame.csv"
    _write_rows(rows, out_csv)
    if errors:
        (output_dir / "errors.log").write_text("\n".join(errors) + "\n", encoding="utf-8")
    if not rows:
        raise RuntimeError(f"No pKa rows were extracted; see {output_dir / 'errors.log'} if present")
    return out_csv
