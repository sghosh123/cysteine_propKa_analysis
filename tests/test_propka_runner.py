from pathlib import Path
import csv
import os
import sys
import pytest

from cysteine_protonation_analysis.propka_runner import run_propka_frames


def test_runner_isolates_frames_and_extracts_requested_residue(tmp_path):
    frames = tmp_path / "frames"
    frames.mkdir()
    (frames / "frame_000000.pdb").write_text("END\n", encoding="utf-8")
    (frames / "frame_000001.pdb").write_text("END\n", encoding="utf-8")
    fake_propka = tmp_path / "fake_propka.py"
    fake_propka.write_text(
        "#!/usr/bin/env python3\n"
        "from pathlib import Path\nimport sys\n"
        "stem = Path(sys.argv[1]).stem\n"
        "Path(stem + '.pka').write_text('CYS 18 P 8.50 0 %\\nCYS 21 P 11.00 0 %\\n')\n",
        encoding="utf-8",
    )
    fake_propka.chmod(0o755)
    output_dir = tmp_path / "run"
    out = run_propka_frames(frames, output_dir, executable=str(fake_propka), resid=18, chain="P")
    assert out.exists()
    with out.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 2
    assert [r["frame_index_zero_based"] for r in rows] == ["0", "1"]
    assert all(r["resid"] == "18" and r["pka"] == "8.5" for r in rows)
    assert (output_dir / "reports/frame_000000/frame_000000.pka").exists()
