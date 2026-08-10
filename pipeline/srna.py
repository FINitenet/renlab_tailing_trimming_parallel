import gzip
import logging
from collections import Counter
from pathlib import Path

import pandas as pd
import pysam
from Bio import SeqIO, bgzf
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .common import run_command, run_parallel


UMI_LINKER = "AACTGTAGGCACCATCAAT"


def _trim(input_file, output_dir, sample, threads, trim_galore, adapter, dry_run):
    command = [
        trim_galore, "--fastqc", "--fastqc_args", f"-t {threads} --nogroup",
        "--basename", sample, "--gzip", "--length", "10", "--trim-n",
        "--suppress_warn", "-j", threads, "-o", output_dir,
    ]
    if adapter:
        command.extend(["-a", adapter, "--max_length", "35", "--consider_already_trimmed", "10"])
    command.append(input_file)
    run_command(command, dry_run=dry_run)


def _extract_umi(trimmed_fastq, output_fasta, sample):
    counts = Counter()
    with pysam.FastqFile(str(trimmed_fastq)) as reads:
        for record in reads:
            if record.sequence.count(UMI_LINKER) != 1:
                continue
            insert, umi = record.sequence.split(UMI_LINKER, 1)
            if 10 <= len(insert) <= 30 and len(umi) == 12 and "N" not in umi:
                counts[(insert, umi)] += 1

    with bgzf.BgzfWriter(str(output_fasta), "wb") as target:
        records = (
            SeqRecord(Seq(insert), id=f"{insert}_{umi}_{count}", description="")
            for (insert, umi), count in counts.items()
        )
        SeqIO.write(records, target, "fasta")


def _map_umi(umi_fasta, umi_dir, sample, config):
    command = [
        config["bowtie"], "-f", "-p", config["threads"], "-m", "50", "-v", "0",
        "--best", "--strata", "-a", "--no-unal", "-x", config["genome_index"],
        umi_fasta, "-S", umi_dir / f"{sample}_aligned.sam", "--al",
        umi_dir / f"{sample}_aligned.fa",
    ]
    if config["dry_run"]:
        run_command(command, dry_run=True)
        return
    result = run_command(command, stdout=-1, stderr=-1)
    (umi_dir / f"{sample}.mapresults.txt").write_bytes(result.stderr)


def _write_length_distribution(fasta_file, output_file):
    lengths = []
    with pysam.FastaFile(str(fasta_file)) as fasta:
        for name in fasta.references:
            lengths.append(fasta.get_reference_length(name))
    counts = pd.Series(lengths, dtype="int64").value_counts().sort_index()
    total = int(counts.sum())
    table = pd.DataFrame({
        "Length": counts.index,
        "Count": counts.values,
        "Percentage": (counts.values / total * 100) if total else [],
    })
    table.to_csv(output_file, index=False, sep="\t")


def _map_reads(trimmed_fastq, mapping_dir, sample, config):
    filtered = mapping_dir / f"{sample}_trimmed_len.fq.gz"
    if config["dry_run"]:
        logging.info("Would retain 18-28 nt reads: %s -> %s", trimmed_fastq, filtered)
    else:
        with gzip.open(trimmed_fastq, "rt") as source, gzip.open(filtered, "wt") as target:
            records = (record for record in SeqIO.parse(source, "fastq") if 18 <= len(record.seq) <= 28)
            SeqIO.write(records, target, "fastq")

    bowtie_command = [
        config["bowtie"], "-p", config["threads"], "-m", "50", "-v", "0",
        "--best", "--strata", "-a", "--no-unal", "-x", config["genome_index"],
        filtered, "-S", mapping_dir / f"{sample}_aligned.sam", "--al",
        mapping_dir / f"{sample}_aligned.fq",
    ]
    if config["dry_run"]:
        run_command(bowtie_command, dry_run=True)
        return
    result = run_command(bowtie_command, stdout=-1, stderr=-1)
    (mapping_dir / f"{sample}.mapresults.txt").write_bytes(result.stderr)


def run_sample(item, config):
    sample, input_file = item
    trim_dir = Path(config["srna_trim_dir"])
    umi_dir = Path(config["srna_umi_dir"])
    mapping_dir = Path(config["srna_mapping_dir"])
    dry_run = config["dry_run"]
    marker = mapping_dir / f"{sample}.mapresults.txt"
    if config["resume"] and marker.is_file():
        logging.info("[%s] sRNA workflow already complete; skipping", sample)
        return
    if not dry_run:
        trim_dir.mkdir(parents=True, exist_ok=True)
        mapping_dir.mkdir(parents=True, exist_ok=True)

    if config["umi_flag"] == 1:
        if not dry_run:
            umi_dir.mkdir(parents=True, exist_ok=True)
        _trim(input_file, umi_dir, sample, config["threads"], config["trim_galore"], None, dry_run)
        umi_trimmed = umi_dir / f"{sample}_trimmed.fq.gz"
        umi_fasta = umi_dir / f"{sample}_umi.fa.gz"
        if dry_run:
            logging.info("Would extract UMI sequences from %s", umi_trimmed)
        else:
            _extract_umi(umi_trimmed, umi_fasta, sample)
        _map_umi(umi_fasta, umi_dir, sample, config)
        if not dry_run:
            _write_length_distribution(
                umi_dir / f"{sample}_aligned.fa", umi_dir / f"{sample}_umi_dist.csv"
            )
        _trim(input_file, trim_dir, sample, config["threads"], config["trim_galore"], UMI_LINKER, dry_run)
    else:
        _trim(input_file, trim_dir, sample, config["threads"], config["trim_galore"], None, dry_run)

    _map_reads(trim_dir / f"{sample}_trimmed.fq.gz", mapping_dir, sample, config)
    logging.info("[%s] sRNA workflow complete", sample)


def run(samples, config, jobs):
    run_parallel(run_sample, samples, config, jobs, "sRNA workflow")
