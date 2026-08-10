#!/usr/bin/env python3
import argparse
import logging
import sys
from pathlib import Path

from pipeline import bubble, length, mapping_summary, preprocess, profile, srna, tailbase
from pipeline.common import (
    discover_fastqs,
    require_bowtie_index,
    require_executable,
    require_file,
)


ROOT = Path(__file__).resolve().parent
ASSETS = ROOT / "assets"
RESOURCES = ROOT / "resources"

DEFAULT_HAIRPIN = (
    "/bios-store1/chenyc/Reference_Source/Arabidopsis_Reference/"
    "ath_hairpin_bowtie_index/hairpin"
)
DEFAULT_TRSNO = (
    "/bios-store1/chenyc/Reference_Source/Arabidopsis_Reference/"
    "ath_trsnoRNA_bowtie_index/ensembl39_araport2016_trsRNA"
)
DEFAULT_GENOME = (
    "/bios-store1/chenyc/Reference_Source/Arabidopsis_Reference/"
    "ath_chr_bowtie_index/Arabidopsis_thaliana.TAIR10.dna.toplevel"
)
DEFAULT_RSCRIPT = "/usr/local/bin/Rscript"
STEP_ORDER = ("preprocess", "profile", "bubble", "srna", "mapping_summary", "length", "tailbase")
STEP_ALIASES = {
    "preprocess": "preprocess", "profile": "profile", "srna": "srna",
    "bubble": "bubble", "bubble_plot": "bubble",
    "srna_workflow": "srna", "mapping": "mapping_summary",
    "mapping_summary": "mapping_summary", "length": "length",
    "length_dist": "length", "tailbase": "tailbase", "tail_base": "tailbase",
}


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Self-contained, parallel plant miRNA tailing/trimming workflow."
    )
    parser.add_argument("-i", "--input", default="1_rawdata", help="raw FASTQ directory")
    parser.add_argument("-o", "--output", default=".", help="project output directory")
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
                        help="miRNA hairpin Bowtie index prefix")
    parser.add_argument("--trsno", default=DEFAULT_TRSNO,
                        help="tRNA/snoRNA Bowtie index prefix")
    parser.add_argument("--genome-index", default=DEFAULT_GENOME,
                        help="genome Bowtie index prefix")
    parser.add_argument("--meta-file", default=str(RESOURCES / "miRNA_start_sequence_length.txt"))
    parser.add_argument("--mechanism-file", default=str(RESOURCES / "ath_miRNA_Mechanism_hairpin.txt"))
    parser.add_argument("--sequence-merge-file", default=str(RESOURCES / "miRNA_sequence_merge.txt"))
    parser.add_argument("--rscript", default=DEFAULT_RSCRIPT,
                        help=f"Rscript executable used by tailbase (default: {DEFAULT_RSCRIPT})")
    parser.add_argument("--bubble-output", default=None,
                        help="multi-page bubble plot PDF; default: <output>/5_GMC_analysis/plot_bubble/tailing_trimming_bubble.pdf")
    parser.add_argument("--bubble-cols", type=int, default=2,
                        help="number of sample panels per row in the bubble plot (default: 2)")
    parser.add_argument("--mapping-tag", default="4_mapping")
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


def _discover_existing_samples(output_dir):
    remapping = output_dir / "3_remapping"
    if remapping.is_dir():
        names = sorted(path.name for path in remapping.iterdir() if path.is_dir())
        if names:
            return [(name, None) for name in names]
    raise FileNotFoundError(
        "Cannot discover samples: provide --filelist, raw input files, or 3_remapping sample directories"
    )


def build_samples(args, steps, output_dir):
    if not set(steps).intersection({"preprocess", "profile", "bubble", "srna", "length", "tailbase"}):
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
    return _discover_existing_samples(output_dir)


def preflight(args, steps, config):
    for key in ("meta_file", "tailbase_script", "mechanism_file", "sequence_merge_file"):
        require_file(config[key], key.replace("_", " "))
    if args.jobs < 1 or args.threads_per_sample < 1 or args.bubble_cols < 1:
        raise ValueError("--jobs, --threads-per-sample, and --bubble-cols must be at least 1")
    if "preprocess" in steps:
        config["trim_galore"] = require_executable("trim_galore")
        config["bowtie"] = require_executable("bowtie")
        require_bowtie_index(args.mir_hairpin)
        require_bowtie_index(args.trsno)
    if "srna" in steps:
        config["trim_galore"] = require_executable("trim_galore")
        config["bowtie"] = require_executable("bowtie")
        require_bowtie_index(args.genome_index)
    if "tailbase" in steps:
        config["rscript"] = require_executable(args.rscript)


def build_config(args):
    output_dir = Path(args.output).expanduser().resolve()
    bubble_output = (
        Path(args.bubble_output).expanduser().resolve() if args.bubble_output
        else output_dir / "5_GMC_analysis" / "plot_bubble" / "tailing_trimming_bubble.pdf"
    )
    mapping_results = (
        Path(args.mapping_results).expanduser().resolve() if args.mapping_results
        else output_dir / f"mapping_results_bowtie_{args.mapping_tag}.csv"
    )
    return {
        "output_dir": output_dir,
        "threads": args.threads_per_sample,
        "adapter": args.adapter,
        "umi_flag": args.umi_flag,
        "mir_hairpin": args.mir_hairpin,
        "trsno": args.trsno,
        "genome_index": args.genome_index,
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


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
    args = parse_args(argv)
    try:
        steps = parse_steps(args.steps)
        config = build_config(args)
        samples = build_samples(args, steps, config["output_dir"])
        preflight(args, steps, config)
        if not args.dry_run:
            config["output_dir"].mkdir(parents=True, exist_ok=True)
        logging.info("Samples (%d): %s", len(samples), ", ".join(item[0] for item in samples))
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
                input_dir = config["output_dir"] / args.mapping_tag
                if args.dry_run:
                    logging.info("Would summarize mapping logs: %s -> %s",
                                 input_dir, config["mapping_results"])
                else:
                    result = mapping_summary.summarize(input_dir, config["mapping_results"])
                    logging.info("Mapping summary contains %d samples", len(result))
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
