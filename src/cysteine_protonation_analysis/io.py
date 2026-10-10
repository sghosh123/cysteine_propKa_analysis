"""Input readers for legacy PROPKA time-series exports and reports."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import numpy as np


@dataclass(frozen=True)
class PKaSeries:
    """A frame-labelled series of predicted pKa values.

    ``frames`` are labels supplied by the source file; they are not assumed to be
    physical time. Most legacy Tcl outputs use one-based frame labels, but users
    should verify each dataset's provenance before converting labels to time.
    """

    frames: np.ndarray
    values: np.ndarray
    source: Path


def read_pka_series(path: str | Path) -> PKaSeries:
    """Read a whitespace-delimited ``frame pKa`` file with strict validation."""
    path = Path(path)
    frames: list[int] = []
    values: list[float] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith(("#", "!")):
                continue
            fields = stripped.split()
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_number}: expected at least two columns (frame pKa)")
            try:
                frame_float = float(fields[0])
                pka_value = float(fields[1])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: non-numeric frame or pKa value") from exc
            if not np.isfinite(frame_float) or not frame_float.is_integer():
                raise ValueError(f"{path}:{line_number}: frame label must be a finite integer")
            if not np.isfinite(pka_value):
                raise ValueError(f"{path}:{line_number}: pKa value must be finite")
            frames.append(int(frame_float))
            values.append(pka_value)
    if not values:
        raise ValueError(f"{path}: no pKa rows found")
    frame_array = np.asarray(frames, dtype=np.int64)
    value_array = np.asarray(values, dtype=np.float64)
    if len(np.unique(frame_array)) != len(frame_array):
        raise ValueError(f"{path}: duplicate frame labels found")
    if np.any(np.diff(frame_array) <= 0):
        raise ValueError(f"{path}: frame labels must be strictly increasing")
    return PKaSeries(frame_array, value_array, path)


_CYS_LINE = re.compile(
    r"^\s*CYS\s+(-?\d+)\s+(\S+)\s+([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)\b"
)


def parse_propka_cysteines(path: str | Path) -> list[dict[str, object]]:
    """Extract primary CYS pKa rows from a PROPKA ``.pka`` text report.

    Interaction-continuation rows are excluded, and repeated residue/chain
    entries in both the detailed and summary report tables are deduplicated.
    """
    path = Path(path)
    rows: list[dict[str, object]] = []
    seen: set[tuple[int, str]] = set()
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            # PROPKA continuation rows repeat the residue and chain but put
            # interaction terms far to the right; limit matching to the leading
            # residue/pKa columns so a later 0.00 contribution is not misread.
            match = _CYS_LINE.match(line[:40])
            if not match:
                continue
            residue_id, chain_id, pka = match.groups()
            key = (int(residue_id), chain_id)
            # Legacy PROPKA reports may repeat the same primary table in the
            # report body and in the final prediction summary. Keep its first
            # occurrence rather than emitting duplicate entries.
            if key in seen:
                continue
            seen.add(key)
            rows.append({
                "source": str(path),
                "line_number": line_number,
                "resname": "CYS",
                "resid": int(residue_id),
                "chain": chain_id,
                "pka": float(pka),
            })
    return rows
