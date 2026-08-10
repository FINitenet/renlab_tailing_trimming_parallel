import logging
from pathlib import Path

from .common import run_command, run_parallel


def run_sample(item, config):
    sample, _ = item
    table_dir = Path(config["tailbase_doc_dir"]) / sample
    plot_dir = Path(config["tailbase_plot_dir"]) / sample
    workspace_dir = Path(config["tailbase_bin_dir"])
    marker = table_dir / f"{sample}_tail_summary.xlsx"
    if config["resume"] and marker.is_file():
        logging.info("[%s] tailbase already complete; skipping", sample)
        return
    if not config["dry_run"]:
        table_dir.mkdir(parents=True, exist_ok=True)
        plot_dir.mkdir(parents=True, exist_ok=True)
        workspace_dir.mkdir(parents=True, exist_ok=True)
    command = [
        config["rscript"], "--vanilla", config["tailbase_script"], sample,
        Path(config["gmc_dir"]) / f"{sample}.5GMC",
        config["mapping_results"], config["meta_file"],
        config["mechanism_file"], config["sequence_merge_file"],
        Path(config["profile_dir"]) / f"{sample}.summary.txt",
        table_dir, plot_dir, workspace_dir,
    ]
    run_command(command, dry_run=config["dry_run"])


def run(samples, config, jobs):
    run_parallel(run_sample, samples, config, jobs, "Tailbase")
