import csv
import logging
import re
from contextlib import ExitStack
from collections import defaultdict
from pathlib import Path

import numpy as np

from .common import run_parallel


def _load_mirnas(meta_file):
    definitions = []
    lookup = {}
    with open(meta_file) as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) < 3:
                continue
            mirna_id, start = row[0], row[2]
            parts = re.split(r"[-.]", mirna_id)
            hairpin = "-".join(parts[:2])
            lookup.setdefault((hairpin, start), mirna_id)
            definitions.append({
                "id": mirna_id,
                "hairpin": hairpin,
                "start": int(start),
                "end": int(start) + int(row[3]) - 1,
                "matrix": [[0] * 11 for _ in range(11)],
            })
    return definitions, lookup


def _parse_query_name(query_name):
    prefix, trim_text = query_name.rsplit("--", 1)
    _, _, sequence = prefix.rsplit("-", 2)
    return sequence, int(trim_text)


def _make_profile_and_5gmc(meta_file, sam_file, profile_file, gmc_file,
                           seqlogo_root, sample):
    definitions, lookup = _load_mirnas(meta_file)
    by_hairpin = defaultdict(list)
    for definition in definitions:
        by_hairpin[definition["hairpin"]].append(definition)

    sample_seqlogo = seqlogo_root / sample
    sample_seqlogo.mkdir(parents=True, exist_ok=True)
    populated = set()
    with ExitStack() as stack, open(sam_file) as source, open(gmc_file, "w") as gmc:
        seqlogo_handles = {}
        gmc.write("RNAME\tID\tStart\tSequence\n")
        for row in csv.reader(source, delimiter="\t"):
            if not row or row[0].startswith("@") or len(row) < 10:
                continue
            reference = row[2]
            position = int(row[3])
            sequence, trim_count = _parse_query_name(row[0])

            mirna_id = lookup.get((reference, row[3]))
            if mirna_id is not None:
                raw_sequence = row[9]
                gmc.write(f"{row[0]}\t{mirna_id}\t{row[3]}\t{raw_sequence}\n")
                if mirna_id not in seqlogo_handles:
                    path = sample_seqlogo / (mirna_id.replace("*", "_star") + ".txt")
                    seqlogo_handles[mirna_id] = stack.enter_context(open(path, "w"))
                seqlogo_handles[mirna_id].write(raw_sequence.replace("T", "U") + "\n")
                populated.add(mirna_id)

            for definition in by_hairpin.get(reference, ()):
                if position != definition["start"]:
                    continue
                end_5gmc = position + len(sequence) - trim_count - 1
                adjusted_trim = trim_count
                if end_5gmc > definition["end"]:
                    end_5gmc = definition["end"]
                    adjusted_trim = position + len(sequence) - definition["end"] - 1
                tailing = definition["end"] - end_5gmc
                if 0 <= adjusted_trim <= 10 and 0 <= tailing <= 10:
                    definition["matrix"][adjusted_trim][tailing] += 1

    for definition in definitions:
        if definition["id"] not in populated:
            path = sample_seqlogo / (definition["id"].replace("*", "_star") + ".txt")
            path.write_text("N" * 30 + "\n")

    with open(profile_file, "w") as target:
        for definition in definitions:
            matrix = definition["matrix"]
            for trim_count in range(10, -1, -1):
                target.write("\t".join(str(matrix[trim_count][tail])
                                       for tail in range(10, -1, -1)) + "\t\n")
            total = sum(sum(row) for row in matrix)
            target.write(f"{definition['id']}\t{total}\t0\n")


def _split_profile(profile_file, output_dir, sample):
    with open(profile_file) as handle:
        matrix = [line.rstrip("\n").split("\t") for line in handle]
    sample_dir = output_dir / sample
    sample_dir.mkdir(parents=True, exist_ok=True)
    for offset in range(0, len(matrix), 12):
        group = matrix[offset:offset + 12]
        if len(group) != 12:
            raise ValueError(f"Incomplete 12-line profile block in {profile_file} at line {offset + 1}")
        tag = group[-1]
        filename = tag[0].split()[0].replace("*", "_star") + ".txt"
        reordered = [tag] + group[:-1]
        with open(sample_dir / filename, "w") as target:
            target.write("\t\n".join("\t".join(row) for row in reordered) + "\t\n")


def _summarize_profile(profile_file, summary_file):
    line_count = sum(1 for _ in open(profile_file))
    if line_count % 12:
        raise ValueError(f"Profile line count is not divisible by 12: {profile_file}")
    result = []
    for index in range(line_count // 12):
        matrix = np.loadtxt(profile_file, skiprows=index * 12, max_rows=11)
        tag = np.loadtxt(profile_file, skiprows=index * 12 + 11, max_rows=1, dtype=str)[0]
        total = np.sum(matrix)
        row = [tag, total, total - np.sum(matrix[:, 10]), total - np.sum(matrix[10, :])]
        row.extend(np.sum(matrix[i, :]) for i in range(9, -1, -1))
        result.append(row)
    np.savetxt(
        summary_file, np.asarray(result), fmt="%s", delimiter="\t",
        header="ID\tSUM\ttrimming\ttailling\ttail_1\ttail_2\ttail_3\ttail_4\t"
               "tail_5\ttail_6\ttail_7\ttail_8\ttail_9\ttail_10",
    )


def run_sample(item, config):
    sample, _ = item
    output_dir = Path(config["output_dir"])
    sam_file = output_dir / "3_remapping" / sample / "merged-alignment-sorted.sam"
    profile_dir = output_dir / "4_163.results"
    summary_dir = output_dir / "5_GMC_analysis"
    profile_file = profile_dir / f"{sample}.txt"
    summary_file = profile_dir / f"{sample}.summary.txt"
    gmc_file = summary_dir / f"{sample}.5GMC"
    if config["resume"] and summary_file.is_file() and gmc_file.is_file():
        logging.info("[%s] profile already complete; skipping", sample)
        return
    if config["dry_run"]:
        logging.info("Would scan %s once and write profile/5GMC outputs for %s", sam_file, sample)
        return
    if not sam_file.is_file():
        raise FileNotFoundError(f"Remapping SAM not found for {sample}: {sam_file}")
    profile_dir.mkdir(parents=True, exist_ok=True)
    summary_dir.mkdir(parents=True, exist_ok=True)
    _make_profile_and_5gmc(config["meta_file"], sam_file, profile_file, gmc_file,
                           summary_dir / "doc_seqlogo", sample)
    _summarize_profile(profile_file, summary_file)
    _split_profile(profile_file, profile_dir, sample)
    logging.info("[%s] profile complete", sample)


def run(samples, config, jobs):
    run_parallel(run_sample, samples, config, jobs, "Profile")
