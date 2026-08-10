import logging
from pathlib import Path

from .common import run_command, run_parallel


def run_sample(item, config):
    sample, _ = item
    output_dir = Path(config["output_dir"])
    marker = output_dir / "5_GMC_analysis" / "doc" / sample / f"{sample}_tail_summary.xlsx"
    if config["resume"] and marker.is_file():
        logging.info("[%s] tailbase already complete; skipping", sample)
        return
    if not config["dry_run"]:
        (output_dir / "5_GMC_analysis" / "doc").mkdir(parents=True, exist_ok=True)
        (output_dir / "5_GMC_analysis" / "plot").mkdir(parents=True, exist_ok=True)
        (output_dir / "5_GMC_analysis" / "bin").mkdir(parents=True, exist_ok=True)
    command = [
        config["rscript"], "--vanilla", config["tailbase_script"], sample,
        output_dir, config["mapping_results"], config["meta_file"],
        config["mechanism_file"], config["sequence_merge_file"],
    ]
    run_command(command, dry_run=config["dry_run"])


def run(samples, config, jobs):
    run_parallel(run_sample, samples, config, jobs, "Tailbase")
