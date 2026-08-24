import csv
import gzip
import logging
import math
import multiprocessing
import os
import re
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/renlab-tailing-matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pysam
from matplotlib.backends.backend_pdf import PdfPages

from .bubble import _plot_sample as _plot_original_sample
from .common import run_command


CANONICAL_TAS_ID = re.compile(r"^TAS(?:1[abc]|2a|3[abc]|4)$", re.IGNORECASE)
COLORS = (
    "#9b4b8f", "#4c9a51", "#4c78a8", "#e6863b", "#b279a2",
    "#54a9a6", "#e45756", "#72b7b2", "#f2cf5b", "#8f6bb3",
)


def _attributes(text):
    return dict(item.split("=", 1) for item in text.split(";") if "=" in item)


def _load_loci(gff_file):
    loci = []
    seen = Counter()
    with open(gff_file) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"gene", "ncRNA_gene"}:
                continue
            attrs = _attributes(fields[8])
            raw_id = attrs.get("ID") or attrs.get("Name") or (
                f"{fields[0]}:{fields[3]}-{fields[4]}"
            )
            biotype = attrs.get("biotype")
            if biotype == "hc-siRNA":
                locus_type = "hc-siRNA"
            elif biotype == "phasi_tasiRNA" and CANONICAL_TAS_ID.fullmatch(raw_id):
                locus_type = "TAS"
            else:
                continue
            seen[raw_id] += 1
            locus_id = raw_id if seen[raw_id] == 1 else f"{raw_id}_{seen[raw_id]}"
            expected = 24 if locus_type == "hc-siRNA" else 21
            loci.append({
                "key": f"L{len(loci) + 1:05d}",
                "id": locus_id,
                "type": locus_type,
                "chrom": fields[0],
                "start": int(fields[3]),
                "end": int(fields[4]),
                "expected": expected,
            })
    if not loci:
        raise ValueError(f"No canonical TAS or hc-siRNA loci found in {gff_file}")
    return loci


def _write_locus_reference(loci, genome_fasta, work_dir, bowtie_build, threads, dry_run):
    fasta_path = work_dir / "TAS_hc_siRNA_loci.fa"
    index_prefix = work_dir / "TAS_hc_siRNA_loci"
    metadata_path = work_dir / "loci.tsv"
    if not dry_run and not fasta_path.is_file():
        with pysam.FastaFile(str(genome_fasta)) as genome, \
                open(fasta_path, "w") as fasta, open(metadata_path, "w") as metadata:
            writer = csv.DictWriter(
                metadata, fieldnames=("key", "id", "type", "chrom", "start", "end", "expected"),
                delimiter="\t",
            )
            writer.writeheader()
            for locus in loci:
                sequence = genome.fetch(locus["chrom"], locus["start"] - 1, locus["end"]).upper()
                fasta.write(f">{locus['key']}\n")
                for offset in range(0, len(sequence), 80):
                    fasta.write(sequence[offset:offset + 80] + "\n")
                writer.writerow(locus)
    index_marker = Path(str(index_prefix) + ".1.ebwt")
    if not index_marker.is_file():
        run_command([bowtie_build, "--threads", threads, fasta_path, index_prefix], dry_run=dry_run)
    return index_prefix


def _open_database(path, rebuild=False):
    if rebuild and path.exists():
        path.unlink()
    connection = sqlite3.connect(path)
    connection.execute("PRAGMA journal_mode=WAL")
    connection.execute("PRAGMA synchronous=NORMAL")
    connection.executescript("""
        CREATE TABLE IF NOT EXISTS sequences (sequence TEXT PRIMARY KEY);
        CREATE TABLE IF NOT EXISTS events (
            sequence TEXT PRIMARY KEY,
            locus_key TEXT NOT NULL,
            strand TEXT NOT NULL,
            anchor INTEGER NOT NULL,
            core_length INTEGER NOT NULL,
            tail_length INTEGER NOT NULL
        );
        CREATE INDEX IF NOT EXISTS events_anchor ON events(locus_key, strand, anchor);
    """)
    return connection


def _count_sequences(samples, mapping_dir, database, work_dir, resume):
    counts_dir = work_dir / "sample_counts"
    counts_dir.mkdir(parents=True, exist_ok=True)
    count_files = [counts_dir / f"{sample}.tsv.gz" for sample, _ in samples]
    unique_fasta = work_dir / "unique_sequences.fa"
    if resume and unique_fasta.is_file() and all(path.is_file() for path in count_files):
        logging.info("Reusing collapsed counts for %d samples", len(samples))
        return unique_fasta, count_files
    for sample_index, (sample, _) in enumerate(samples):
        count_file = counts_dir / f"{sample}.tsv.gz"
        if resume and count_file.is_file():
            with gzip.open(count_file, "rt") as handle:
                database.executemany(
                    "INSERT OR IGNORE INTO sequences(sequence) VALUES (?)",
                    ((line.split("\t", 1)[0],) for line in handle),
                )
            database.commit()
            continue
        fastq = Path(mapping_dir) / f"{sample}_trimmed_len.fq.gz"
        if not fastq.is_file():
            raise FileNotFoundError(f"Filtered 18-28 nt FASTQ not found for {sample}: {fastq}")
        counts = Counter()
        with gzip.open(fastq, "rt") as handle:
            while True:
                name = handle.readline()
                if not name:
                    break
                sequence = handle.readline().strip().upper().replace("U", "T")
                plus = handle.readline()
                quality = handle.readline()
                if not plus or not quality:
                    raise ValueError(f"Truncated FASTQ record in {fastq}")
                if 18 <= len(sequence) <= 28 and set(sequence) <= set("ACGT"):
                    counts[sequence] += 1
        with gzip.open(count_file, "wt") as target:
            for sequence, count in sorted(counts.items()):
                target.write(f"{sequence}\t{count}\n")
        database.executemany(
            "INSERT OR IGNORE INTO sequences(sequence) VALUES (?)",
            ((sequence,) for sequence in counts),
        )
        database.commit()
        logging.info("[%s] collapsed %d distinct 18-28 nt sequences", sample, len(counts))
    if not unique_fasta.is_file() or not resume:
        with open(unique_fasta, "w") as target:
            for (sequence,) in database.execute("SELECT sequence FROM sequences ORDER BY sequence"):
                target.write(f">{sequence}\n{sequence}\n")
    return unique_fasta, count_files


def _trim_unmapped(source_path, target_path, trim_length):
    with open(source_path) as source, open(target_path, "w") as target:
        name = None
        chunks = []
        for line in source:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    sequence = "".join(chunks)
                    if len(sequence) > 12:
                        target.write(f">{name}\n{sequence[:-1]}\n")
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            sequence = "".join(chunks)
            if len(sequence) > 12:
                target.write(f">{name}\n{sequence[:-1]}\n")


def _import_sam_events(sam_path, tail_length, database):
    rows = []
    with pysam.AlignmentFile(str(sam_path), "r") as alignments:
        for record in alignments.fetch(until_eof=True):
            if record.is_unmapped:
                continue
            core_length = record.query_alignment_length
            anchor = (
                record.reference_start + 1 if not record.is_reverse
                else record.reference_end
            )
            rows.append((
                record.query_name, record.reference_name,
                "-" if record.is_reverse else "+", anchor, core_length, tail_length,
            ))
            if len(rows) >= 100000:
                database.executemany(
                    "INSERT OR IGNORE INTO events VALUES (?, ?, ?, ?, ?, ?)", rows
                )
                rows.clear()
    if rows:
        database.executemany("INSERT OR IGNORE INTO events VALUES (?, ?, ?, ?, ?, ?)", rows)
    database.commit()


def _map_unique_sequences(unique_fasta, index_prefix, database, work_dir,
                          bowtie, threads, resume, dry_run):
    current = unique_fasta
    for tail_length in range(11):
        sam_path = work_dir / f"mapped_tail_{tail_length}.sam"
        unmapped = work_dir / f"unmapped_tail_{tail_length}.fa"
        marker = work_dir / f"mapped_tail_{tail_length}.done"
        if not (resume and marker.is_file()):
            run_command([
                bowtie, "-f", "-p", threads, "-v", "0", "--best",
                "-k", "1", "--no-unal", "--un", unmapped, "-S",
                index_prefix, current, sam_path,
            ], dry_run=dry_run)
            if dry_run:
                current = unmapped
                continue
            _import_sam_events(sam_path, tail_length, database)
            marker.touch()
        if tail_length < 10:
            next_input = work_dir / f"trimmed_tail_{tail_length + 1}.fa"
            if not (resume and next_input.is_file()):
                _trim_unmapped(unmapped, next_input, tail_length + 1)
            current = next_input


def _load_event_map(database):
    return {
        sequence: (locus_key, strand, anchor, core_length, tail_length)
        for sequence, locus_key, strand, anchor, core_length, tail_length
        in database.execute("SELECT sequence, locus_key, strand, anchor, core_length, tail_length FROM events")
    }


def _iter_counts(count_file):
    with gzip.open(count_file, "rt") as handle:
        for line in handle:
            sequence, count = line.rstrip("\n").split("\t")
            yield sequence, int(count)


def _canonical_lengths(event_map, count_files, loci_by_key):
    pooled = defaultdict(Counter)
    for count_file in count_files:
        for sequence, count in _iter_counts(count_file):
            event = event_map.get(sequence)
            if event is None or event[4] != 0:
                continue
            locus_key, strand, anchor, core_length, _ = event
            pooled[(locus_key, strand, anchor)][core_length] += count
    canonical = {}
    for anchor, length_counts in pooled.items():
        expected = loci_by_key[anchor[0]]["expected"]
        canonical[anchor] = max(
            length_counts.items(), key=lambda value: (value[1], -abs(value[0] - expected))
        )[0]
    return canonical


def _aggregate(event_map, count_files, samples, loci_by_key, canonical):
    matrices = defaultdict(lambda: np.zeros((11, 11), dtype=float))
    totals = defaultdict(float)
    exact_totals = defaultdict(float)
    for sample_index, count_file in enumerate(count_files):
        for sequence, count in _iter_counts(count_file):
            event = event_map.get(sequence)
            if event is None:
                continue
            locus_key, strand, anchor, core_length, tail_length = event
            key = (locus_key, strand, anchor)
            reference_length = canonical.get(key)
            if reference_length is None:
                continue
            trimming = max(0, reference_length - core_length)
            if trimming > 10 or tail_length > 10:
                continue
            matrices[(locus_key, sample_index)][trimming, tail_length] += count
            totals[(locus_key, sample_index)] += count
            if tail_length == 0:
                exact_totals[(locus_key, sample_index)] += count
    return matrices, totals, exact_totals


def _selected_loci(loci, totals, exact_totals, sample_count, hc_min_abundance):
    selected = []
    for locus in loci:
        per_sample = [totals.get((locus["key"], index), 0) for index in range(sample_count)]
        exact_per_sample = [
            exact_totals.get((locus["key"], index), 0) for index in range(sample_count)
        ]
        maximum_exact = max(exact_per_sample, default=0)
        if locus["type"] == "TAS" or maximum_exact >= hc_min_abundance:
            selected.append((locus, maximum_exact, sum(per_sample)))
    return selected


def _write_tables(output_dir, selected, matrices, totals, samples, canonical):
    with open(output_dir / "selected_loci.tsv", "w") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("locus_id", "type", "chrom", "start", "end", "max_sample_exact_abundance", "pooled_supported_abundance"))
        for locus, maximum, pooled in selected:
            writer.writerow((locus["id"], locus["type"], locus["chrom"], locus["start"], locus["end"], maximum, pooled))
    with open(output_dir / "bubble_matrix_long.tsv", "w") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("locus_id", "type", "sample", "trimming", "tailing", "count", "percentage"))
        for locus, _, _ in selected:
            for sample_index, (sample, _) in enumerate(samples):
                matrix = matrices.get((locus["key"], sample_index), np.zeros((11, 11)))
                total = totals.get((locus["key"], sample_index), 0)
                for trimming, tailing in zip(*np.nonzero(matrix)):
                    count = matrix[trimming, tailing]
                    writer.writerow((locus["id"], locus["type"], sample, trimming, tailing,
                                     count, count / total * 100 if total else 0))
    with open(output_dir / "local_anchor_reference_lengths.tsv", "w") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("locus_key", "strand", "five_prime_anchor", "canonical_length"))
        for key, length in sorted(canonical.items()):
            writer.writerow((*key, length))


def _plot_panel(axis, matrix, sample, color):
    # bubble._plot_sample reads the historical profile matrix in descending
    # trimming/tailing order. Reverse both dimensions so locus matrices use
    # exactly the same plotting code and appearance as the original workflow.
    display_matrix = matrix[::-1, ::-1]
    _plot_original_sample(axis, display_matrix, sample, int(matrix.sum()), color)


def _plot_pdf(path, selected, locus_type, matrices, sample_group, columns):
    loci = [item[0] for item in selected if item[0]["type"] == locus_type]
    rows = math.ceil(len(sample_group) / columns)
    with PdfPages(path) as pdf:
        for index, locus in enumerate(loci, start=1):
            figure, axes = plt.subplots(
                rows, columns, squeeze=False,
                figsize=(5.0 * columns, 4.8 * rows),
            )
            for panel_index, (sample_index, sample) in enumerate(sample_group):
                matrix = matrices.get((locus["key"], sample_index), np.zeros((11, 11)))
                _plot_panel(
                    axes.flat[panel_index], matrix, sample,
                    COLORS[panel_index % len(COLORS)],
                )
            for axis in axes.flat[len(sample_group):]:
                axis.set_visible(False)
            figure.suptitle(
                f"{locus['id']} ({locus_type}; {locus['chrom']}:{locus['start']}-{locus['end']})",
                fontsize=12, fontweight="bold",
            )
            figure.tight_layout(rect=(0, 0, 1, 0.975))
            pdf.savefig(figure)
            plt.close(figure)
            if index == 1 or index % 25 == 0 or index == len(loci):
                logging.info("%s bubble pages: %d/%d", locus_type, index, len(loci))


def _sample_groups(samples, metadata_path, group_columns):
    index_by_raw = {sample: index for index, (sample, _) in enumerate(samples)}
    if metadata_path is None:
        return {"All": [(index, sample) for index, (sample, _) in enumerate(samples)]}
    if not group_columns:
        raise ValueError("At least one locus bubble group column is required")
    groups = defaultdict(list)
    with open(metadata_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample_raw", *group_columns}
        if not required <= set(reader.fieldnames or ()):
            missing_columns = sorted(required - set(reader.fieldnames or ()))
            raise ValueError(
                f"Metadata is missing locus bubble columns {missing_columns}: {metadata_path}"
            )
        for row in reader:
            raw = row["sample_raw"]
            if raw in index_by_raw:
                values = []
                for column in group_columns:
                    value = str(row[column]).strip()
                    if column.lower() in {"replicate", "rep", "repeat"}:
                        value = f"R{value}"
                    values.append(value)
                group_name = "_".join(values)
                groups[group_name].append((index_by_raw[raw], raw))
    missing = set(index_by_raw) - {sample for group in groups.values() for _, sample in group}
    if missing:
        raise ValueError("Samples missing from locus bubble metadata: " + ", ".join(sorted(missing)))
    return dict(groups)


def _plot_tas_group(output_dir, group_name, sample_group, selected, matrices, config):
    tas_pdf = output_dir / f"{group_name}_TAS_tailing_trimming_bubble.pdf"
    _plot_pdf(
        tas_pdf, selected, "TAS", matrices, sample_group,
        config["locus_bubble_cols"],
    )


def _plot_groups(output_dir, sample_groups, selected, matrices, config):
    # Each PDF is independent. Forking lets custom metadata produce separate
    # replicate/group files concurrently without changing the original panel
    # drawing routine or copying the large matrix dictionary at process start.
    context = multiprocessing.get_context("fork")

    def run_processes(tasks):
        failed = []
        max_parallel = min(48, os.cpu_count() or 1)
        for offset in range(0, len(tasks), max_parallel):
            processes = []
            for name, target, args in tasks[offset:offset + max_parallel]:
                process = context.Process(target=target, args=args, name=name)
                process.start()
                processes.append(process)
            for process in processes:
                process.join()
                if process.exitcode != 0:
                    failed.append(f"{process.name} (exit {process.exitcode})")
        if failed:
            raise RuntimeError("Locus bubble plotting failed: " + ", ".join(failed))

    groups = list(sample_groups.items())
    tas_tasks = []
    for group_name, sample_group in groups:
        tas_tasks.append((
            f"TAS-{group_name}", _plot_tas_group,
            (output_dir, group_name, sample_group, selected, matrices, config),
        ))
    run_processes(tas_tasks)

    hc_selected = [item for item in selected if item[0]["type"] == "hc-siRNA"]
    if not hc_selected:
        raise ValueError("No hc-siRNA loci meet the configured abundance threshold")
    chunk_count = min(8, len(hc_selected))
    chunk_size = math.ceil(len(hc_selected) / chunk_count)
    parts_dir = output_dir / ".hc_pdf_parts"
    parts_dir.mkdir(parents=True, exist_ok=True)
    hc_tasks = []
    parts_by_group = defaultdict(list)
    for group_name, sample_group in groups:
        for chunk_index in range(chunk_count):
            start = chunk_index * chunk_size
            chunk = hc_selected[start:start + chunk_size]
            if not chunk:
                continue
            part_path = parts_dir / f"{group_name}.part_{chunk_index + 1:02d}.pdf"
            parts_by_group[group_name].append(part_path)
            hc_tasks.append((
                f"hc-{group_name}-{chunk_index + 1}", _plot_pdf,
                (part_path, chunk, "hc-siRNA", matrices, sample_group,
                 config["locus_bubble_cols"]),
            ))
    run_processes(hc_tasks)

    for group_name, _ in groups:
        hc_pdf = output_dir / (
            f"{group_name}_hc_siRNA_abundance_ge_"
            f"{config['locus_bubble_hc_min_abundance']}_tailing_trimming_bubble.pdf"
        )
        merged_pdf = parts_dir / f"{group_name}.merged.pdf"
        if merged_pdf.exists():
            merged_pdf.unlink()
        run_command([config["pdfunite"], *parts_by_group[group_name], merged_pdf])
        os.replace(merged_pdf, hc_pdf)
        for part_path in parts_by_group[group_name]:
            part_path.unlink()
    parts_dir.rmdir()


def run(samples, config):
    if not samples:
        raise ValueError("No samples selected for locus tailing/trimming bubbles")
    output_dir = Path(config["locus_bubble_output"])
    work_dir = Path(config["locus_bubble_work_dir"])
    sample_groups = _sample_groups(
        samples, config["locus_bubble_metadata"], config["locus_bubble_group_columns"]
    )
    completion_marker = output_dir / "canonical_TAS_per_locus_hc.complete"
    expected_pdfs = []
    for tissue in sample_groups:
        expected_pdfs.extend((
            output_dir / f"{tissue}_TAS_tailing_trimming_bubble.pdf",
            output_dir / (
                f"{tissue}_hc_siRNA_abundance_ge_"
                f"{config['locus_bubble_hc_min_abundance']}_tailing_trimming_bubble.pdf"
            ),
        ))
    if (config["resume"] and completion_marker.is_file()
            and all(path.is_file() for path in expected_pdfs)):
        logging.info("Locus bubble plots already complete; skipping: %s", output_dir)
        return
    if config["dry_run"]:
        logging.info("Would build TAS and abundance-filtered hc-siRNA locus bubble plots in %s", output_dir)
    else:
        output_dir.mkdir(parents=True, exist_ok=True)
        work_dir.mkdir(parents=True, exist_ok=True)

    # Keep locus-definition-dependent mappings separate from the legacy run,
    # which included PHAS21/PHAS24 loci under the TAS label. Collapsed sample
    # counts stay in the parent work directory and can still be reused.
    mapping_work_dir = work_dir / "canonical_TAS_v2"
    if not config["dry_run"]:
        mapping_work_dir.mkdir(parents=True, exist_ok=True)

    loci = _load_loci(config["rnatype_annotation"])
    loci_by_key = {locus["key"]: locus for locus in loci}
    index_prefix = _write_locus_reference(
        loci, config["genome_fasta"], mapping_work_dir, config["bowtie_build"],
        config["threads"], config["dry_run"],
    )
    if config["dry_run"]:
        return
    database_path = mapping_work_dir / "locus_bubble.sqlite"
    database = _open_database(database_path)
    try:
        unique_fasta, count_files = _count_sequences(
            samples, config["srna_mapping_dir"], database, work_dir, config["resume"]
        )
        _map_unique_sequences(
            unique_fasta, index_prefix, database, mapping_work_dir, config["bowtie"],
            config["threads"], config["resume"], config["dry_run"],
        )
        event_map = _load_event_map(database)
        logging.info("Loaded %d locus-mapped distinct sequences", len(event_map))
        canonical = _canonical_lengths(event_map, count_files, loci_by_key)
        matrices, totals, exact_totals = _aggregate(
            event_map, count_files, samples, loci_by_key, canonical
        )
        del event_map
    finally:
        database.close()

    selected = _selected_loci(
        loci, totals, exact_totals, len(samples), config["locus_bubble_hc_min_abundance"]
    )
    _write_tables(output_dir, selected, matrices, totals, samples, canonical)
    _plot_groups(output_dir, sample_groups, selected, matrices, config)
    logging.info(
        "Locus bubbles complete: %d TAS loci and %d hc-siRNA loci",
        sum(item[0]["type"] == "TAS" for item in selected),
        sum(item[0]["type"] == "hc-siRNA" for item in selected),
    )
    completion_marker.touch()
