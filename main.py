#!/usr/bin/env python3
import argparse
import logging
import os
import sys
from pathlib import Path

from pipeline import bubble, length, mapping_summary, preprocess, profile, rnatype, srna, tailbase
from pipeline.common import (
    discover_fastqs,
    require_bowtie_index,
    require_executable,
    require_file,
)


ROOT = Path(__file__).resolve().parent
ASSETS = ROOT / "assets"
RESOURCES = ROOT / "resources"
VERSION = "0.1.0"

DEFAULT_HAIRPIN = os.environ.get("MIR3END_HAIRPIN_INDEX")
DEFAULT_TRSNO = os.environ.get("MIR3END_TRSNO_INDEX")
DEFAULT_GENOME = os.environ.get("MIR3END_GENOME_INDEX")
DEFAULT_GENOME_FASTA = os.environ.get("MIR3END_GENOME_FASTA")
DEFAULT_RNATYPE_ANNOTATION = os.environ.get("MIR3END_RNATYPE_ANNOTATION")
DEFAULT_SHORTSTACK = os.environ.get("MIR3END_SHORTSTACK", "ShortStack")
DEFAULT_FEATURECOUNTS = os.environ.get("MIR3END_FEATURECOUNTS", "featureCounts")
DEFAULT_RSCRIPT = os.environ.get("MIR3END_RSCRIPT", "Rscript")
STEP_ORDER = (
    "preprocess", "profile", "bubble", "srna", "mapping_summary", "rnatype", "length", "tailbase"
)
STEP_ALIASES = {
    "preprocess": "preprocess", "profile": "profile", "srna": "srna",
    "bubble": "bubble", "bubble_plot": "bubble",
    "srna_workflow": "srna", "mapping": "mapping_summary",
    "mapping_summary": "mapping_summary", "length": "length",
    "rnatype": "rnatype", "rna_type": "rnatype", "rna-type": "rnatype",
    "length_dist": "length", "tailbase": "tailbase", "tail_base": "tailbase",
}


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Self-contained, parallel plant miRNA tailing/trimming workflow."
    )
    parser.add_argument("--version", action="version", version=f"miR3End {VERSION}")
    parser.add_argument("-i", "--input", default="1_rawdata", help="raw FASTQ directory")
    parser.add_argument("-o", "--output", default=".", help="project output directory")
    parser.add_argument("--layout", choices=("auto", "organized", "legacy"), default="auto",
                        help="output layout: auto preserves existing legacy projects; new projects use organized")
    parser.add_argument("-j", "--jobs", type=int, default=1,
                        help="maximum concurrent samples (default: 1)")
    parser.add_argument("--threads-per-sample", type=int, default=8,
                        help="threads used by each external tool process (default: 8)")
    parser.add_argument("--suffix", default=".fastq.gz",
                        help="input filename suffix removed to form sample names")
    parser.add_argument("--adapter", default="AGATCGGAAGAG")
    parser.add_argument("--umi-flag", type=int, choices=(1, 2), default=2,
                        help="1: UMI library; 2: standard small-RNA library")
    parser.add_argument("--mir-hairpin", default=DEFAULT_HAIRPIN,
                        help="miRNA hairpin Bowtie index prefix (or MIR3END_HAIRPIN_INDEX)")
    parser.add_argument("--trsno", default=DEFAULT_TRSNO,
                        help="tRNA/snoRNA Bowtie index prefix (or MIR3END_TRSNO_INDEX)")
    parser.add_argument("--genome-index", default=DEFAULT_GENOME,
                        help="genome Bowtie index prefix (or MIR3END_GENOME_INDEX)")
    parser.add_argument("--genome-fasta", default=DEFAULT_GENOME_FASTA,
                        help="genome FASTA used by ShortStack (or MIR3END_GENOME_FASTA)")
    parser.add_argument("--rnatype-annotation", default=DEFAULT_RNATYPE_ANNOTATION,
                        help="GFF3 biotype annotation (or MIR3END_RNATYPE_ANNOTATION)")
    parser.add_argument("--rnatype-priority", default=None,
                        help="optional RNA-type priority file, one category per line")
    parser.add_argument("--rnatype-stranded", type=int, choices=(0, 1, 2), default=0,
                        help="featureCounts strand mode: 0 unstranded, 1 forward, 2 reverse")
    parser.add_argument("--shortstack", default=DEFAULT_SHORTSTACK,
                        help="ShortStack executable")
    parser.add_argument("--featurecounts", default=DEFAULT_FEATURECOUNTS,
                        help=f"featureCounts executable (default: {DEFAULT_FEATURECOUNTS})")
    parser.add_argument("--meta-file", default=str(RESOURCES / "miRNA_start_sequence_length.txt"))
    parser.add_argument("--mechanism-file", default=str(RESOURCES / "ath_miRNA_Mechanism_hairpin.txt"))
    parser.add_argument("--sequence-merge-file", default=str(RESOURCES / "miRNA_sequence_merge.txt"))
    parser.add_argument("--rscript", default=DEFAULT_RSCRIPT,
                        help=f"Rscript executable used by tailbase (default: {DEFAULT_RSCRIPT})")
    parser.add_argument("--bubble-output", default=None,
                        help="multi-page bubble plot PDF; default depends on --layout")
    parser.add_argument("--bubble-cols", type=int, default=2,
                        help="number of sample panels per row in the bubble plot (default: 2)")
    parser.add_argument("--mapping-tag", default="4_mapping",
                        help="legacy mapping directory/output tag (default: 4_mapping)")
    parser.add_argument("--mapping-results", default=None,
                        help="mapping summary table; default: <output>/mapping_results_bowtie_<tag>.csv")
    parser.add_argument("--filelist", default=None,
                        help="optional sample list; otherwise samples are discovered automatically")
    parser.add_argument("--steps", default="all",
                        help="comma-separated steps or all: " + ",".join(STEP_ORDER))
    parser.add_argument("--resume", action="store_true",
                        help="skip a sample step when its final output already exists")
    parser.add_argument("--dry-run", action="store_true", help="validate and print without executing")
    return parser.parse_args(argv)


def parse_steps(raw):
    if raw.strip() == "all":
        return list(STEP_ORDER)
    steps = []
    for value in raw.split(","):
        key = value.strip()
        if not key:
            continue
        if key not in STEP_ALIASES:
            raise ValueError(f"Unknown step '{key}'. Allowed: all, {', '.join(STEP_ORDER)}")
        step = STEP_ALIASES[key]
        if step not in steps:
            steps.append(step)
    if not steps:
        raise ValueError("No steps selected")
    return steps


def _read_filelist(path):
    with open(path) as handle:
        names = [line.strip() for line in handle if line.strip()]
    if not names:
        raise ValueError(f"No samples in filelist: {path}")
    return [(name.replace("-", "_"), None) for name in names]


def _discover_existing_samples(remapping):
    if remapping.is_dir():
        names = sorted(path.name for path in remapping.iterdir() if path.is_dir())
        if names:
            return [(name, None) for name in names]
    raise FileNotFoundError(
        f"Cannot discover samples: provide --filelist, raw input files, or sample directories in {remapping}"
    )


def build_samples(args, steps, config):
    if not set(steps).intersection(
        {"preprocess", "profile", "bubble", "srna", "rnatype", "length", "tailbase"}
    ):
        return []
    if args.filelist:
        return _read_filelist(Path(args.filelist).expanduser().resolve())
    input_dir = Path(args.input).expanduser().resolve()
    if input_dir.is_dir():
        try:
            return discover_fastqs(input_dir, args.suffix)
        except FileNotFoundError:
            if "preprocess" in steps or "srna" in steps:
                raise
    return _discover_existing_samples(config["remapping_dir"])


def preflight(args, steps, config):
    if set(steps).intersection({"profile", "bubble", "length", "tailbase"}):
        require_file(config["meta_file"], "meta file")
    if "tailbase" in steps:
        for key in ("tailbase_script", "mechanism_file", "sequence_merge_file"):
            require_file(config[key], key.replace("_", " "))
    if args.jobs < 1 or args.threads_per_sample < 1 or args.bubble_cols < 1:
        raise ValueError("--jobs, --threads-per-sample, and --bubble-cols must be at least 1")
    if "preprocess" in steps:
        config["trim_galore"] = require_executable("trim_galore")
        config["bowtie"] = require_executable("bowtie")
        if not args.mir_hairpin or not args.trsno:
            raise ValueError(
                "--mir-hairpin and --trsno are required for preprocess "
                "(or set MIR3END_HAIRPIN_INDEX and MIR3END_TRSNO_INDEX)"
            )
        require_bowtie_index(args.mir_hairpin)
        require_bowtie_index(args.trsno)
    if "srna" in steps:
        config["trim_galore"] = require_executable("trim_galore")
        config["bowtie"] = require_executable("bowtie")
        if not args.genome_index:
            raise ValueError(
                "--genome-index is required for srna (or set MIR3END_GENOME_INDEX)"
            )
        require_bowtie_index(args.genome_index)
    if "rnatype" in steps:
        config["shortstack"] = require_executable(args.shortstack)
        config["featurecounts"] = require_executable(args.featurecounts)
        if not args.genome_fasta or not args.rnatype_annotation:
            raise ValueError(
                "--genome-fasta and --rnatype-annotation are required for rnatype "
                "(or set MIR3END_GENOME_FASTA and MIR3END_RNATYPE_ANNOTATION)"
            )
        require_file(config["genome_fasta"], "ShortStack genome FASTA")
        require_file(config["rnatype_annotation"], "RNA-type GFF3 annotation")
        if config["rnatype_priority"] is not None:
            require_file(config["rnatype_priority"], "RNA-type priority file")
    if "tailbase" in steps:
        config["rscript"] = require_executable(args.rscript)


def build_config(args):
    output_dir = Path(args.output).expanduser().resolve()
    layout = args.layout
    if layout == "auto":
        legacy_markers = (output_dir / "3_remapping", output_dir / "4_163.results")
        organized_markers = (output_dir / "work", output_dir / "results")
        if any(path.exists() for path in organized_markers):
            layout = "organized"
        elif any(path.exists() for path in legacy_markers):
            layout = "legacy"
        else:
            layout = "organized"

    if layout == "organized":
        tt_work = output_dir / "work" / "tailing_trimming"
        srna_work = output_dir / "work" / "srna"
        tt_results = output_dir / "results" / "tailing_trimming"
        srna_results = output_dir / "results" / "srna"
        paths = {
            "format_dir": tt_work / "trimmed_fasta",
            "remapping_dir": tt_work / "remapping",
            "profile_dir": tt_results / "profiles",
            "gmc_dir": tt_results / "5gmc",
            "seqlogo_dir": tt_results / "sequence_logo",
            "length_dir": tt_results / "length",
            "tailbase_doc_dir": tt_results / "tailbase" / "tables",
            "tailbase_plot_dir": tt_results / "tailbase" / "plots",
            "tailbase_bin_dir": tt_results / "tailbase" / "workspace",
            "srna_trim_dir": srna_work / "trimmed",
            "srna_umi_dir": srna_work / "umi",
            "srna_mapping_dir": srna_work / "genome_mapping",
            "rnatype_shortstack_dir": srna_work / "rnatype" / "shortstack",
            "rnatype_annotation_dir": srna_work / "rnatype" / "featurecounts",
            "rnatype_results_dir": srna_results / "rnatype",
        }
        default_bubble = tt_results / "bubble" / "tailing_trimming_bubble.pdf"
        default_mapping_results = srna_results / "mapping_summary.tsv"
    else:
        summary_dir = output_dir / "5_GMC_analysis"
        paths = {
            "format_dir": output_dir / "2_formmat_fasta",
            "remapping_dir": output_dir / "3_remapping",
            "profile_dir": output_dir / "4_163.results",
            "gmc_dir": summary_dir,
            "seqlogo_dir": summary_dir / "doc_seqlogo",
            "length_dir": output_dir,
            "tailbase_doc_dir": summary_dir / "doc",
            "tailbase_plot_dir": summary_dir / "plot",
            "tailbase_bin_dir": summary_dir / "bin",
            "srna_trim_dir": output_dir / "2_trim_adapter",
            "srna_umi_dir": output_dir / "3_extract_umi",
            "srna_mapping_dir": output_dir / "4_mapping",
            "rnatype_shortstack_dir": output_dir / "6_RNA_type_analysis" / "work" / "shortstack",
            "rnatype_annotation_dir": output_dir / "6_RNA_type_analysis" / "work" / "featurecounts",
            "rnatype_results_dir": output_dir / "6_RNA_type_analysis" / "results",
        }
        default_bubble = summary_dir / "plot_bubble" / "tailing_trimming_bubble.pdf"
        default_mapping_results = output_dir / f"mapping_results_bowtie_{args.mapping_tag}.csv"

    bubble_output = (
        Path(args.bubble_output).expanduser().resolve() if args.bubble_output
        else default_bubble
    )
    mapping_results = (
        Path(args.mapping_results).expanduser().resolve() if args.mapping_results
        else default_mapping_results
    )
    config = {
        "output_dir": output_dir,
        "layout": layout,
        "threads": args.threads_per_sample,
        "adapter": args.adapter,
        "umi_flag": args.umi_flag,
        "mir_hairpin": args.mir_hairpin,
        "trsno": args.trsno,
        "genome_index": args.genome_index,
        "genome_fasta": (
            Path(args.genome_fasta).expanduser().resolve() if args.genome_fasta else None
        ),
        "rnatype_annotation": (
            Path(args.rnatype_annotation).expanduser().resolve()
            if args.rnatype_annotation else None
        ),
        "rnatype_priority": (
            Path(args.rnatype_priority).expanduser().resolve() if args.rnatype_priority else None
        ),
        "rnatype_stranded": args.rnatype_stranded,
        "meta_file": Path(args.meta_file).expanduser().resolve(),
        "mechanism_file": Path(args.mechanism_file).expanduser().resolve(),
        "sequence_merge_file": Path(args.sequence_merge_file).expanduser().resolve(),
        "tailbase_script": ASSETS / "tail_base_summary.R",
        "bubble_output": bubble_output,
        "bubble_cols": args.bubble_cols,
        "mapping_results": mapping_results,
        "resume": args.resume,
        "dry_run": args.dry_run,
    }
    config.update(paths)
    return config


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
    args = parse_args(argv)
    try:
        steps = parse_steps(args.steps)
        config = build_config(args)
        samples = build_samples(args, steps, config)
        preflight(args, steps, config)
        if not args.dry_run:
            config["output_dir"].mkdir(parents=True, exist_ok=True)
        logging.info("Samples (%d): %s", len(samples), ", ".join(item[0] for item in samples))
        logging.info("Output layout: %s", config["layout"])
        logging.info("Steps: %s; jobs=%d; threads/sample=%d",
                     ", ".join(steps), args.jobs, args.threads_per_sample)

        for step in steps:
            logging.info("Starting step: %s", step)
            if step == "preprocess":
                preprocess.run(samples, config, args.jobs)
            elif step == "profile":
                profile.run(samples, config, args.jobs)
            elif step == "bubble":
                bubble.run(samples, config)
            elif step == "srna":
                srna.run(samples, config, args.jobs)
            elif step == "mapping_summary":
                input_dir = config["srna_mapping_dir"]
                if args.dry_run:
                    logging.info("Would summarize mapping logs: %s -> %s",
                                 input_dir, config["mapping_results"])
                else:
                    result = mapping_summary.summarize(input_dir, config["mapping_results"])
                    logging.info("Mapping summary contains %d samples", len(result))
            elif step == "rnatype":
                rnatype.run(samples, config, args.jobs)
            elif step == "length":
                length.run(samples, config, args.jobs)
            elif step == "tailbase":
                if not args.dry_run:
                    require_file(config["mapping_results"], "mapping results")
                tailbase.run(samples, config, args.jobs)
            logging.info("Finished step: %s", step)
        logging.info("Pipeline completed successfully")
        return 0
    except Exception as error:
        logging.exception("Pipeline failed: %s", error)
        return 1


if __name__ == "__main__":
    sys.exit(main())
