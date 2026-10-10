from pathlib import Path
from cysteine_protonation_analysis.provenance import audit_terminal_values


def test_terminal_value_audit_flags_mismatch_and_match(tmp_path):
    series_root = tmp_path / "series"
    reports_root = tmp_path / "reports"
    (series_root / "system").mkdir(parents=True)
    (reports_root / "system").mkdir(parents=True)
    (series_root / "system/pKa_Cys18P_test.dat").write_text("1 9.0\n2 9.5\n", encoding="utf-8")
    (reports_root / "system/out.pka").write_text("CYS 18 P 9.5 0 %\n", encoding="utf-8")
    (series_root / "system/pKa_Cys21P_test.dat").write_text("1 10.0\n2 11.5\n", encoding="utf-8")
    (reports_root / "system/out.pka").write_text("CYS 18 P 9.5 0 %\nCYS 21 P 11.0 0 %\n", encoding="utf-8")
    rows = audit_terminal_values(series_root, reports_root)
    by_name = {Path(row["series"]).name: row for row in rows}
    assert by_name["pKa_Cys18P_test.dat"]["status"] == "match"
    assert by_name["pKa_Cys21P_test.dat"]["status"] == "mismatch"
    assert by_name["pKa_Cys21P_test.dat"]["difference_last_minus_report"] == 0.5
