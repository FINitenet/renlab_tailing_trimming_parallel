import gzip
import logging
import os
import shutil
import uuid
from collections import Counter, OrderedDict
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/renlab-tailing-matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import pysam
from Bio import SeqIO

from .common import run_command, run_parallel


DEFAULT_PRIORITY = OrderedDict((
    ("miRNA_primary_transcript", 1),
    ("phasi_tasiRNA", 2),
    ("phasi/tasiRNA", 2),
    ("hc-siRNA", 3),
    ("transposable_element", 4),
    ("lncRNA", 5),
    ("protein_coding", 6),
    ("snoRNA", 7),
    ("snRNA", 8),
    ("rRNA", 9),
    ("tRNA", 10),
    ("otherRNA", 11),
    ("unassigned", 11),
))

DISPLAY_ORDER = (
    "unassigned", "snoRNA", "snRNA", "rRNA", "tRNA", "lncRNA",
    "protein_coding", "transposable_element", "hc-siRNA",
    "phasi_tasiRNA", "miRNA_primary_transcript",
)

COLORS = {
    "unassigned": "#a6cee3",
    "snoRNA": "#33a02c",
    "snRNA": "#fb9a99",
    "rRNA": "#1f78b4",
    "tRNA": "#b2df8a",
    "lncRNA": "#e31a1c",
    "protein_coding": "#fdbf6f",
    "transposable_element": "#ff7f00",
    "hc-siRNA": "#cab2d6",
    "phasi_tasiRNA": "#6a3d9a",
    "miRNA_primary_transcript": "#ffff99",
}


def _priority(path):
    if path is None:
        return dict(DEFAULT_PRIORITY)
    result = {}
    with open(path) as handle:
        for rank, line in enumerate(handle, 1):
            name = line.strip()
            if name and not name.startswith("#") and name not in result:
                result[name] = rank
    if not result:
        raise ValueError(f"RNA-type priority file contains no categories: {path}")
    result.setdefault("otherRNA", len(result) + 1)
    result.setdefault("unassigned", result["otherRNA"])
    return result


def _ensure_filtered_reads(sample, config):
    filtered = Path(config["srna_mapping_dir"]) / f"{sample}_trimmed_len.fq.gz"
    if filtered.is_file():
        return filtered
    trimmed = Path(config["srna_trim_dir"]) / f"{sample}_trimmed.fq.gz"
    if not trimmed.is_file():
        raise FileNotFoundError(
            f"Length-filtered reads not found for {sample}: {filtered}. "
            "Run --steps srna,rnatype or provide the existing sRNA work directory."
        )
    logging.info("[%s] retaining 18-28 nt reads for RNA-type analysis", sample)
    filtered.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(trimmed, "rt") as source, gzip.open(filtered, "wt") as target:
        reads = (record for record in SeqIO.parse(source, "fastq") if 18 <= len(record.seq) <= 28)
        SeqIO.write(reads, target, "fastq")
    return filtered


def _run_shortstack(sample, reads, config):
    shortstack_dir = Path(config["rnatype_shortstack_dir"])
    final_bam = shortstack_dir / f"{sample}.bam"
    if final_bam.is_file():
        return final_bam
    shortstack_dir.mkdir(parents=True, exist_ok=True)
    temporary = shortstack_dir / (
        f".{sample}.shortstack.{os.getpid()}.{uuid.uuid4().hex[:8]}"
    )
    command = [
        config["shortstack"],
        "--genomefile", config["genome_fasta"],
        "--outdir", temporary,
        "--align_only", "--nohp", "--keep_quals",
        "--mmap", "u", "--bowtie_m", "1000", "--ranmax", "50",
        "--mismatches", "0", "--bowtie_cores", config["threads"],
        "--readfile", reads,
    ]
    run_command(command)
    bam_files = list(temporary.glob("*.bam"))
    if len(bam_files) != 1:
        raise RuntimeError(
            f"ShortStack produced {len(bam_files)} BAM files for {sample}; expected exactly one"
        )
    shutil.move(str(bam_files[0]), final_bam)
    log_file = temporary / "Log.txt"
    if log_file.is_file():
        shutil.move(str(log_file), shortstack_dir / f"{sample}.log")
    shutil.rmtree(temporary)
    return final_bam


def _run_featurecounts(sample, bam_file, config):
    annotation_dir = Path(config["rnatype_annotation_dir"])
    annotation_dir.mkdir(parents=True, exist_ok=True)
    annotation = _prepare_annotation_saf(config["rnatype_annotation"], annotation_dir)
    tagged_bam = annotation_dir / f"{sample}.featureCounts.bam"
    if tagged_bam.is_file():
        return tagged_bam
    count_file = annotation_dir / f"{sample}.featureCounts.txt"
    command = [
        config["featurecounts"],
        "-a", annotation,
        "-F", "SAF",
        "-o", count_file,
        "-T", config["threads"],
        "-s", config["rnatype_stranded"],
        "-O", "-M", "-R", "BAM", "--largestOverlap", "--fraction",
        bam_file,
    ]
    log_file = annotation_dir / f"{sample}.featureCounts.log"
    with open(log_file, "wb") as log_handle:
        run_command(command, stdout=log_handle, stderr=log_handle)
    generated_candidates = (
        annotation_dir / f"{bam_file.name}.featureCounts.bam",
        Path(str(bam_file) + ".featureCounts.bam"),
    )
    generated = next((path for path in generated_candidates if path.is_file()), None)
    if generated is None:
        expected = ", ".join(str(path) for path in generated_candidates)
        raise FileNotFoundError(f"featureCounts tagged BAM not found for {sample}; expected {expected}")
    shutil.move(str(generated), tagged_bam)
    return tagged_bam


def _prepare_annotation_saf(annotation_file, output_dir):
    """Convert gene-level GFF3 records to SAF for featureCounts.

    featureCounts accepts only one value for ``-t``.  The Arabidopsis
    annotation represents RNA loci as both ``gene`` and ``ncRNA_gene``, so a
    comma-separated ``-t gene,ncRNA_gene`` silently selects no records.  SAF
    lets both record types be retained while using the GFF3 ``biotype`` as the
    feature identifier consumed by the downstream classifier.
    """
    annotation_file = Path(annotation_file)
    output = Path(output_dir) / f"{annotation_file.stem}.gene_ncRNA_gene.by_biotype.saf"
    if output.is_file() and output.stat().st_size > 0:
        return output

    temporary = output.with_name(f".{output.name}.{os.getpid()}.{uuid.uuid4().hex[:8]}")
    selected = 0
    with open(annotation_file) as source, open(temporary, "w") as target:
        target.write("GeneID\tChr\tStart\tEnd\tStrand\n")
        for line in source:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"gene", "ncRNA_gene"}:
                continue
            attributes = {}
            for item in fields[8].split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    attributes[key] = value
            biotype = attributes.get("biotype")
            if not biotype:
                continue
            target.write(
                f"{biotype}\t{fields[0]}\t{fields[3]}\t{fields[4]}\t{fields[6]}\n"
            )
            selected += 1
    if not selected:
        temporary.unlink(missing_ok=True)
        raise ValueError(
            f"No gene or ncRNA_gene records with biotype found in annotation: {annotation_file}"
        )
    os.replace(temporary, output)
    return output


def _choose_type(record, priorities):
    status = record.get_tag("XS") if record.has_tag("XS") else None
    if status == "Assigned" and record.has_tag("XT"):
        features = record.get_tag("XT").split(",")
        known = [feature for feature in features if feature in priorities]
        if known:
            selected = min(known, key=priorities.get)
            if selected == "phasi/tasiRNA":
                return "phasi_tasiRNA"
            if selected == "otherRNA":
                return "unassigned"
            return selected
        return "unassigned"
    if status == "Unassigned_NoFeatures":
        return "unassigned"
    return None


def _classify_bam(tagged_bam, priorities):
    counts = Counter()
    with pysam.AlignmentFile(str(tagged_bam), "rb") as alignments:
        for record in alignments.fetch(until_eof=True):
            if record.is_secondary or record.is_supplementary:
                continue
            length = record.query_alignment_length
            if length is None:
                length = record.query_length
            if length is None:
                continue
            rna_type = _choose_type(record, priorities)
            if rna_type is None:
                continue
            counts[(int(length), rna_type)] += 1
    return counts


def _write_sample_table(path, counts):
    path.parent.mkdir(parents=True, exist_ok=True)
    total = sum(counts.values())
    rows = []
    for (length, rna_type), count in sorted(counts.items()):
        rows.append({
            "Length": length,
            "RNA_type": rna_type,
            "Count": count,
            "Percentage": count / total * 100 if total else 0,
        })
    pd.DataFrame(rows, columns=("Length", "RNA_type", "Count", "Percentage")).to_csv(
        path, sep="\t", index=False, float_format="%.6f"
    )


def run_sample(item, config):
    sample, _ = item
    table = Path(config["rnatype_results_dir"]) / "tables" / f"{sample}_rnatype_by_length.tsv"
    if config["resume"] and table.is_file():
        logging.info("[%s] RNA-type classification already complete; skipping", sample)
        return
    if config["dry_run"]:
        reads = Path(config["srna_mapping_dir"]) / f"{sample}_trimmed_len.fq.gz"
        logging.info(
            "Would classify RNA types: %s -> ShortStack -> featureCounts -> %s", reads, table
        )
        return
    reads = _ensure_filtered_reads(sample, config)
    bam_file = _run_shortstack(sample, reads, config)
    tagged_bam = _run_featurecounts(sample, bam_file, config)
    counts = _classify_bam(tagged_bam, _priority(config["rnatype_priority"]))
    _write_sample_table(table, counts)
    logging.info("[%s] RNA-type classification complete (%d mapped reads)", sample, sum(counts.values()))


def _sample_summary(samples, tables_dir):
    count_rows = []
    length_rows = []
    for sample, _ in samples:
        path = tables_dir / f"{sample}_rnatype_by_length.tsv"
        if not path.is_file():
            raise FileNotFoundError(f"RNA-type table not found for {sample}: {path}")
        table = pd.read_csv(path, sep="\t")
        grouped = table.groupby("RNA_type", as_index=False)["Count"].sum()
        totals = grouped["Count"].sum()
        for row in grouped.itertuples(index=False):
            count_rows.append({
                "Sample": sample,
                "RNA_type": row.RNA_type,
                "Count": int(row.Count),
                "Percentage": row.Count / totals * 100 if totals else 0,
            })
        sample_table = table.copy()
        sample_table.insert(0, "Sample", sample)
        length_rows.append(sample_table)
    return pd.DataFrame(count_rows), pd.concat(length_rows, ignore_index=True)


def _plot_overall(summary, output_base):
    matrix = summary.pivot(index="Sample", columns="RNA_type", values="Percentage").fillna(0)
    categories = [name for name in DISPLAY_ORDER if name in matrix.columns]
    categories.extend(name for name in matrix.columns if name not in categories)
    matrix = matrix[categories]
    width = max(10, len(matrix) * 0.42)
    figure, axis = plt.subplots(figsize=(width, 6))
    bottom = pd.Series(0.0, index=matrix.index)
    # Draw in reverse so the visible stack from top to bottom follows the
    # legend order: unassigned -> ... -> miRNA_primary_transcript.
    for category in reversed(categories):
        values = matrix[category]
        axis.bar(
            matrix.index, values, bottom=bottom, label=category,
            color=COLORS.get(category, "#808080"), width=0.82,
        )
        bottom = bottom + values
    axis.set_ylabel("Percentage (%)")
    axis.set_xlabel("")
    axis.set_ylim(0, 100)
    axis.tick_params(axis="x", labelrotation=90)
    handles, labels = axis.get_legend_handles_labels()
    axis.legend(
        handles[::-1], labels[::-1], title="RNA type",
        bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False,
    )
    axis.grid(axis="y", color="#d0d0d0", linewidth=0.5, alpha=0.8)
    axis.set_axisbelow(True)
    figure.tight_layout()
    figure.savefig(output_base.with_suffix(".pdf"), bbox_inches="tight")
    figure.savefig(output_base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(figure)


def _plot_by_length(table, sample, output_file):
    matrix = table.pivot_table(
        index="Length", columns="RNA_type", values="Percentage", aggfunc="sum", fill_value=0
    )
    categories = [name for name in DISPLAY_ORDER if name in matrix.columns]
    categories.extend(name for name in matrix.columns if name not in categories)
    matrix = matrix[categories]
    figure, axis = plt.subplots(figsize=(9, 5))
    bottom = pd.Series(0.0, index=matrix.index)
    # Draw in reverse so the visible stack from top to bottom follows the
    # legend order: unassigned -> ... -> miRNA_primary_transcript.
    for category in reversed(categories):
        values = matrix[category]
        axis.bar(
            matrix.index, values, bottom=bottom, label=category,
            color=COLORS.get(category, "#808080"), width=0.82,
        )
        bottom = bottom + values
    axis.set_xlabel("Read length (nt)")
    axis.set_ylabel("Percentage (%)")
    axis.set_xticks(matrix.index)
    axis.set_title(sample, fontweight="bold")
    handles, labels = axis.get_legend_handles_labels()
    axis.legend(
        handles[::-1], labels[::-1], title="RNA type",
        bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False,
    )
    axis.grid(axis="y", color="#d0d0d0", linewidth=0.5, alpha=0.8)
    axis.set_axisbelow(True)
    figure.tight_layout()
    figure.savefig(output_file, bbox_inches="tight")
    plt.close(figure)


def summarize(samples, config):
    results_dir = Path(config["rnatype_results_dir"])
    tables_dir = results_dir / "tables"
    plots_dir = results_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    summary, by_length = _sample_summary(samples, tables_dir)
    summary.to_csv(results_dir / "rnatype_summary.tsv", sep="\t", index=False, float_format="%.6f")
    with pd.ExcelWriter(results_dir / "rnatype_summary.xlsx", engine="openpyxl") as writer:
        summary.to_excel(writer, sheet_name="RNA type summary", index=False)
        by_length.to_excel(writer, sheet_name="RNA type by length", index=False)
    _plot_overall(summary, plots_dir / "rnatype_composition")
    for sample, table in by_length.groupby("Sample", sort=False):
        _plot_by_length(table, sample, plots_dir / f"{sample}_rnatype_by_length.pdf")
    logging.info("RNA-type summary complete for %d samples", len(samples))


def run(samples, config, jobs):
    run_parallel(run_sample, samples, config, jobs, "RNA-type classification")
    if config["dry_run"]:
        logging.info("Would aggregate RNA-type tables and plots for %d samples", len(samples))
    else:
        summarize(samples, config)
