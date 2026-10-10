from pathlib import Path
import numpy as np
import pytest

from cysteine_protonation_analysis.io import read_pka_series, parse_propka_cysteines
from cysteine_protonation_analysis.metrics import summarize_series, henderson_hasselbalch_deprotonated_fraction


def test_read_series_and_summary(tmp_path):
    path = tmp_path / "values.dat"
    path.write_text("# frame pKa\n1 6.0\n2 7.0\n3 8.0\n", encoding="utf-8")
    series = read_pka_series(path)
    assert series.frames.tolist() == [1, 2, 3]
    assert series.values.tolist() == [6.0, 7.0, 8.0]
    out = summarize_series(series, ph=7.0)
    assert out["n_frames"] == 3
    assert out["median_pka"] == 7.0
    assert out["fraction_frames_with_pka_below_pH"] == pytest.approx(1 / 3)
    assert out["mean_HH_thiolate_fraction_model_only"] == pytest.approx(0.5)


def test_reject_duplicate_or_unordered_frame_labels(tmp_path):
    path = tmp_path / "bad.dat"
    path.write_text("1 7\n1 8\n", encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate"):
        read_pka_series(path)
    path.write_text("2 7\n1 8\n", encoding="utf-8")
    with pytest.raises(ValueError, match="strictly increasing"):
        read_pka_series(path)


def test_reject_nonfinite_values_and_bad_columns(tmp_path):
    path = tmp_path / "bad.dat"
    path.write_text("1 nan\n", encoding="utf-8")
    with pytest.raises(ValueError, match="finite"):
        read_pka_series(path)
    path.write_text("1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="two columns"):
        read_pka_series(path)


def test_henderson_hasselbalch_equal_at_ph_equals_pka():
    fractions = henderson_hasselbalch_deprotonated_fraction(np.array([7.0]), 7.0)
    assert fractions[0] == pytest.approx(0.5)


def test_parse_propka_primary_cys_rows_only(tmp_path):
    path = tmp_path / "example.pka"
    path.write_text(
        "CYS  18 P   9.47   0 %  0.47 204\n"
        "CYS  18 P                                            0.00 XXX 0 X\n"
        "CYS  21 P  12.25  57 %  3.27 442\n"
        "CYS  18 P   9.47      9.00  # summary duplicate\n",
        encoding="utf-8",
    )
    rows = parse_propka_cysteines(path)
    assert [(r["resid"], r["chain"], r["pka"]) for r in rows] == [(18, "P", 9.47), (21, "P", 12.25)]
