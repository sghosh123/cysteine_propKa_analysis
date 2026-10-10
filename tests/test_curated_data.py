from pathlib import Path
import csv
from collections import Counter

from cysteine_protonation_analysis.io import read_pka_series
from cysteine_protonation_analysis.provenance import audit_terminal_values

ROOT = Path(__file__).resolve().parents[1]


def test_all_curated_pka_series_match_manifest_row_counts():
    with (ROOT / "data/manifest.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 44
    for row in rows:
        series = read_pka_series(ROOT / row["path"])
        assert len(series.values) == int(row["n_rows"])


def test_terminal_audit_preserves_four_known_provenance_mismatches():
    rows = audit_terminal_values(ROOT / "data/raw/propka_series", ROOT / "data/raw/propka_reports")
    counts = Counter(row["status"] for row in rows)
    assert counts["match"] == 40
    assert counts["mismatch"] == 4
    mismatches = {row["series"] for row in rows if row["status"] == "mismatch"}
    assert mismatches == {
        "yeast/y-complex_2ggp/pKa_Cys13C_yeast_MBD_complex.dat",
        "yeast/y-complex_2ggp/pKa_Cys15M_yeast_MBD_complex.dat",
        "yeast/y-complex_2ggp/pKa_Cys16C_yeast_MBD_complex.dat",
        "yeast/y-complex_2ggp/pKa_Cys18M_yeast_MBD_complex.dat",
    }
