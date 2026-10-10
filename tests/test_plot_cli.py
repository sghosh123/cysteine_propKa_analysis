from pathlib import Path
import csv
import json
import subprocess
import sys

from cysteine_protonation_analysis.cli import main
from cysteine_protonation_analysis.io import read_pka_series
from cysteine_protonation_analysis.plotting import plot_series


def test_plot_writes_both_outputs(tmp_path):
    source = tmp_path / "series.dat"
    source.write_text("1 7.1\n2 8.2\n3 9.0\n", encoding="utf-8")
    outputs = plot_series(read_pka_series(source), tmp_path / "figs" / "cys", bins=3)
    assert all(path.exists() and path.stat().st_size > 0 for path in outputs)


def test_cli_summarize_json_and_csv(tmp_path, capsys):
    source = tmp_path / "series.dat"
    source.write_text("1 7.1\n2 8.2\n3 9.0\n", encoding="utf-8")
    assert main(["summarize", str(source), "--ph", "7.0"]) == 0
    data = json.loads(capsys.readouterr().out)
    assert data["n_frames"] == 3
    out = tmp_path / "summary.csv"
    assert main(["summarize", str(source), "--output", str(out)]) == 0
    with out.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 1 and rows[0]["n_frames"] == "3"


def test_cli_parse_report(tmp_path):
    source = tmp_path / "x.pka"
    source.write_text("CYS 18 P 9.00 0 %\n", encoding="utf-8")
    out = tmp_path / "cys.csv"
    assert main(["parse-report", str(source), "--output", str(out)]) == 0
    assert "9.0" in out.read_text(encoding="utf-8")


def test_cli_compare_manifest(tmp_path):
    d1 = tmp_path / "one.dat"
    d2 = tmp_path / "two.dat"
    d1.write_text("1 7\n2 8\n", encoding="utf-8")
    d2.write_text("1 9\n2 10\n", encoding="utf-8")
    manifest = tmp_path / "manifest.csv"
    manifest.write_text(f"path\n{d1}\n{d2}\n", encoding="utf-8")
    output = tmp_path / "summary.csv"
    assert main(["compare", "--manifest", str(manifest), "--output", str(output)]) == 0
    assert len(list(csv.DictReader(output.open()))) == 2
