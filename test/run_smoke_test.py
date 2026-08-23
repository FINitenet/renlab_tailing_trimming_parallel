#!/usr/bin/env python3
"""Run a small, dependency-light end-to-end test of the profile stage."""

from __future__ import annotations

import csv
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "test" / "data"


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="mir3end-smoke-") as temp_name:
        project = Path(temp_name)
        remapping = project / "work" / "tailing_trimming" / "remapping" / "demo"
        remapping.mkdir(parents=True)
        shutil.copy2(
            DATA / "merged-alignment-sorted.sam",
            remapping / "merged-alignment-sorted.sam",
        )

        command = [
            sys.executable,
            str(ROOT / "main.py"),
            "--output", str(project),
            "--layout", "organized",
            "--steps", "profile",
            "--filelist", str(DATA / "samples.txt"),
            "--meta-file", str(DATA / "miRNA_start_sequence_length.txt"),
        ]
        subprocess.run(command, check=True, cwd=ROOT)

        summary = (
            project / "results" / "tailing_trimming" / "profiles" / "demo.summary.txt"
        )
        gmc = project / "results" / "tailing_trimming" / "5gmc" / "demo.5GMC"
        split_profile = (
            project / "results" / "tailing_trimming" / "profiles" / "demo"
            / "ath-miR156a-5p.txt"
        )
        for required in (summary, gmc, split_profile):
            if not required.is_file() or required.stat().st_size == 0:
                raise AssertionError(f"Missing or empty output: {required}")

        with summary.open() as handle:
            rows = list(csv.reader(handle, delimiter="\t"))
        data_rows = [row for row in rows if row and not row[0].startswith("#")]
        if len(data_rows) != 1 or data_rows[0][0] != "ath-miR156a-5p":
            raise AssertionError(f"Unexpected profile summary: {data_rows}")
        if float(data_rows[0][1]) != 3.0:
            raise AssertionError(f"Expected three assigned reads: {data_rows[0]}")

    print("miR3End smoke test passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
