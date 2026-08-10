import re
from pathlib import Path

import pandas as pd


def _first_integer(line):
    match = re.search(r"\d+", line)
    if not match:
        raise ValueError(f"No integer found in mapping log line: {line.rstrip()}")
    return int(match.group())


def summarize(input_dir, output_file):
    rows = []
    for log_file in sorted(Path(input_dir).glob("*.mapresults.txt")):
        values = {}
        for line in log_file.read_text().splitlines():
            if "reads processed" in line:
                values["total"] = _first_integer(line)
            elif "reads with at least one alignment" in line:
                values["mapped"] = _first_integer(line)
            elif "Reported" in line:
                values["reported"] = _first_integer(line)
        missing = {"total", "mapped", "reported"} - values.keys()
        if missing:
            raise ValueError(f"Missing {sorted(missing)} in mapping log: {log_file}")
        rate = round(values["mapped"] / values["total"] * 100, 2) if values["total"] else 0
        rows.append({
            "Sample": log_file.name[:-len(".mapresults.txt")],
            "Total_reads": values["total"],
            "Mapped_reads": values["mapped"],
            "Mapped_rate(%)": rate,
            "Reported": values["reported"],
        })
    if not rows:
        raise FileNotFoundError(f"No '*.mapresults.txt' files found in {input_dir}")
    result = pd.DataFrame(rows)
    result.to_csv(output_file, index=False, sep="\t")
    return result
