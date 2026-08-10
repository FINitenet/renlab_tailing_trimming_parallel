import logging
from collections import Counter, OrderedDict
from pathlib import Path

import pysam
from openpyxl import Workbook

from .common import run_parallel


def _build_end_matrix(meta_file):
    result = {}
    with open(meta_file) as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            parts = fields[0].split("-")
            hairpin = "-".join(parts[:2])
            result[f"{hairpin}-{fields[2]}"] = int(fields[2]) + len(fields[1]) - 1
    return result


def run_sample(item, config):
    sample, _ = item
    output_dir = Path(config["output_dir"])
    sam_file = output_dir / "3_remapping" / sample / "merged-alignment-sorted.sam"
    output_file = output_dir / f"{sample}_len_dist.xlsx"
    if config["resume"] and output_file.is_file():
        logging.info("[%s] length distribution already complete; skipping", sample)
        return
    if config["dry_run"]:
        logging.info("Would calculate length distribution: %s -> %s", sam_file, output_file)
        return
    if not sam_file.is_file():
        raise FileNotFoundError(f"Remapping SAM not found for {sample}: {sam_file}")

    ends = _build_end_matrix(config["meta_file"])
    read_lengths = Counter()
    gmc_lengths = Counter()
    seen = set()
    total = 0
    with pysam.AlignmentFile(str(sam_file), "r") as alignments:
        for record in alignments:
            mirna_id = f"{record.reference_name}-{record.reference_start + 1}"
            name_parts = record.query_name.split("-")
            query_id = "-".join(name_parts[:2])
            query_cut = int(record.query_name.split("--")[1])
            raw_length = len(name_parts[2])
            end_5gmc = record.reference_start + raw_length - query_cut
            length_5gmc = raw_length - query_cut
            if mirna_id in ends and query_id not in seen:
                if query_cut == 0 and end_5gmc > ends[mirna_id]:
                    length_5gmc = raw_length
                read_lengths[raw_length] += 1
                gmc_lengths[length_5gmc] += 1
                total += 1
            seen.add(query_id)

    workbook = Workbook()
    total_sheet = workbook.active
    total_sheet.title = "total reads mapped to miRNA"
    total_sheet.append(("total reads mapped to miRNA:", total))
    read_sheet = workbook.create_sheet("distribution of whole reads")
    gmc_sheet = workbook.create_sheet("distribution of 5GMC")
    for sheet, counts in ((read_sheet, read_lengths), (gmc_sheet, gmc_lengths)):
        sheet.append(("Length", "Count", "Percentage"))
        for length, count in OrderedDict(sorted(counts.items())).items():
            percentage = count / total if total else 0
            sheet.append((length, count, f"{percentage:.4f}"))
    workbook.save(output_file)
    logging.info("[%s] length distribution complete", sample)


def run(samples, config, jobs):
    run_parallel(run_sample, samples, config, jobs, "Length analysis")
