import gzip
import logging
from pathlib import Path

import pysam
from Bio import SeqIO

from .common import run_command, run_parallel


def _rewrite_fastq_as_fasta(trimmed_fastq, output_fasta, sample):
    with gzip.open(trimmed_fastq, "rt") as source, open(output_fasta, "w") as target:
        def records():
            for index, record in enumerate(SeqIO.parse(source, "fastq")):
                record.id = f"{sample}-{index}-score:37"
                record.description = ""
                yield record
        SeqIO.write(records(), target, "fasta")


def _reformat_fasta(input_fasta, output_fasta):
    with open(input_fasta) as source, open(output_fasta, "w") as target:
        def records():
            for record in SeqIO.parse(source, "fasta"):
                record.id = record.id.replace("-score:37", f"-{record.seq}--0")
                record.description = ""
                yield record
        SeqIO.write(records(), target, "fasta")


def _run_bowtie(command, log_handle, dry_run):
    run_command(command, dry_run=dry_run, stdout=log_handle, stderr=log_handle)


def run_sample(item, config):
    sample, input_file = item
    format_dir = Path(config["format_dir"])
    remapping_root = Path(config["remapping_dir"])
    remapping_dir = remapping_root / sample
    threads = config["threads"]
    dry_run = config["dry_run"]

    output_fasta = format_dir / f"{sample}_trimmed.fasta"
    trimmed_fastq = format_dir / f"{sample}_trimmed.fq.gz"
    merged_sam = remapping_dir / "merged-alignment-sorted.sam"
    if config["resume"] and merged_sam.is_file():
        logging.info("[%s] preprocess already complete; skipping", sample)
        return

    if not dry_run:
        format_dir.mkdir(parents=True, exist_ok=True)
        remapping_dir.mkdir(parents=True, exist_ok=True)

    trim_command = [
        config["trim_galore"], "--fastqc", "--fastqc_args", f"-t {threads} --nogroup",
        "--gzip", "-q", "20", "--length", "12", "--max_length", "30",
        "--trim-n", "--stringency", "3", "--basename", sample, "--no_report_file",
        "-a", config["adapter"],
        "-j", threads, "-o", format_dir, input_file,
    ]
    log_path = format_dir / f"{sample}.log.txt"
    if dry_run:
        run_command(trim_command, dry_run=True)
        logging.info("Would convert %s to %s", trimmed_fastq, output_fasta)
    else:
        with open(log_path, "a") as log_handle:
            run_command(trim_command, stderr=log_handle)
        _rewrite_fastq_as_fasta(trimmed_fastq, output_fasta, sample)
        _reformat_fasta(output_fasta, remapping_dir / "unmapped-to-genome.fasta")

    if dry_run:
        log_handle = None
    else:
        log_handle = open(remapping_dir / "remapping.log", "w")
    try:
        unmapped_genome = remapping_dir / "unmapped-to-genome.fasta"
        clean_fasta = remapping_dir / "clean.fasta"
        _run_bowtie([
            config["bowtie"], "-p", threads, "-v", "0", "-S", "-a", "-f",
            "--un", clean_fasta, config["trsno"], unmapped_genome,
            remapping_dir / "map2trsnoRNA.sam",
        ], log_handle, dry_run)
        _run_bowtie([
            config["bowtie"], "-p", threads, "-v", "0", "-S", "-a", "-f",
            "--un", remapping_dir / "unmapped-0.fasta", config["mir_hairpin"],
            clean_fasta, remapping_dir / "mapped-0.sam",
        ], log_handle, dry_run)

        for trim_length in range(1, 11):
            previous = remapping_dir / f"unmapped-{trim_length - 1}.fasta"
            trimmed = remapping_dir / f"trimmed-{trim_length}.fasta"
            if not dry_run:
                with open(previous) as source, open(trimmed, "w") as target:
                    def records():
                        for record in SeqIO.parse(source, "fasta"):
                            record.seq = record.seq[:-1]
                            if len(record.seq) >= 12:
                                record.id = record.id.split("--")[0] + f"--{trim_length}"
                                record.description = ""
                                yield record
                    SeqIO.write(records(), target, "fasta")
            _run_bowtie([
                config["bowtie"], "-p", threads, "-v", "0", "-S", "-a", "--norc",
                "--no-unal", "-f", "--un",
                remapping_dir / f"unmapped-{trim_length}.fasta", config["mir_hairpin"], trimmed,
                remapping_dir / f"mapped-{trim_length}.sam",
            ], log_handle, dry_run)
    finally:
        if log_handle is not None:
            log_handle.close()

    if dry_run:
        logging.info("Would normalize, sort, index, and merge SAM/BAM files for %s", sample)
        return

    sorted_bams = []
    for trim_length in range(11):
        sam_path = remapping_dir / f"mapped-{trim_length}.sam"
        bam_path = remapping_dir / f"mapped-{trim_length}.bam"
        sorted_bam = remapping_dir / f"mapped-{trim_length}.sorted.bam"
        with pysam.AlignmentFile(str(sam_path), "r", threads=threads) as source:
            with pysam.AlignmentFile(str(bam_path), "wb", template=source, threads=threads) as target:
                for record in source:
                    sequence = record.query_name.split("-")[-3]
                    record.query_sequence = sequence
                    record.template_length = len(sequence)
                    record.cigartuples = [(0, len(sequence))]
                    trim_count = int(record.query_name.split("--")[1])
                    record.set_tag("TM", sequence[-trim_count:] if trim_count else "")
                    target.write(record)
        pysam.sort("-@", str(threads), "-o", str(sorted_bam), str(bam_path))
        pysam.index(str(sorted_bam))
        sam_path.unlink()
        sorted_bams.append(sorted_bam)

    pysam.merge("-@", str(threads), "-f", "-O", "SAM", str(merged_sam),
                *(str(path) for path in sorted_bams))
    logging.info("[%s] preprocess/remapping complete", sample)


def run(samples, config, jobs):
    run_parallel(run_sample, samples, config, jobs, "Preprocess")
