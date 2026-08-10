import csv
import logging
import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/renlab-tailing-matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages


COLORS = (
    "#9b4b8f", "#4c9a51", "#4c78a8", "#e6863b", "#b279a2",
    "#54a9a6", "#e45756", "#72b7b2", "#f2cf5b", "#8f6bb3",
)


def _load_mirna_ids(meta_file):
    with open(meta_file) as handle:
        return [row[0] for row in csv.reader(handle, delimiter="\t") if row]


def _read_matrix(path):
    with open(path) as handle:
        rows = [line.rstrip("\n").split("\t") for line in handle]
    if len(rows) < 12:
        raise ValueError(f"Incomplete trimming/tailing matrix: {path}")
    try:
        total = int(float(rows[0][1]))
        matrix = np.asarray(
            [[float(value) for value in row[:11]] for row in rows[1:12]],
            dtype=float,
        )
    except (IndexError, ValueError) as error:
        raise ValueError(f"Invalid trimming/tailing matrix: {path}") from error
    if matrix.shape != (11, 11):
        raise ValueError(f"Expected an 11 x 11 matrix in {path}; found {matrix.shape}")
    return total, matrix


def _plot_sample(axis, matrix, sample, total, color):
    coordinates = np.arange(10, -1, -1)
    x_values, y_values = np.meshgrid(coordinates, coordinates)
    matrix_sum = matrix.sum()
    sizes = matrix.ravel() / matrix_sum * 1500 if matrix_sum else np.zeros(121)

    axis.scatter(
        x_values.ravel(), y_values.ravel(), s=sizes, color=color,
        edgecolors="none", zorder=3,
    )
    axis.set_xlim(11, -1)
    axis.set_ylim(-1, 11)
    axis.set_xticks(range(11))
    axis.set_yticks(range(11))
    axis.grid(True, color="#b8b8b8", linewidth=0.5, zorder=0)
    axis.set_aspect("equal")
    axis.tick_params(length=0, labelsize=8)
    axis.set_xlabel("Length of trimming", fontsize=9, fontweight="bold")
    axis.set_ylabel("Length of tailing", fontsize=9, fontweight="bold")
    axis.set_title(f"{sample}\ntotal reads: {total}", fontsize=9, fontweight="bold")


def run(samples, config):
    output_file = Path(config["bubble_output"])
    sample_names = [sample for sample, _ in samples]
    if not sample_names:
        raise ValueError("No samples selected for bubble plotting")
    if config["resume"] and output_file.is_file():
        logging.info("Bubble plot already complete; skipping: %s", output_file)
        return
    if config["dry_run"]:
        logging.info(
            "Would generate a multi-page trimming/tailing bubble plot for %d samples: %s",
            len(sample_names), output_file,
        )
        return

    profile_root = Path(config["profile_dir"])
    for sample in sample_names:
        sample_dir = profile_root / sample
        if not sample_dir.is_dir():
            raise FileNotFoundError(f"Profile matrix directory not found for {sample}: {sample_dir}")

    mirna_ids = _load_mirna_ids(config["meta_file"])
    if not mirna_ids:
        raise ValueError(f"No miRNA definitions found in {config['meta_file']}")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    columns = min(config["bubble_cols"], len(sample_names))
    rows = math.ceil(len(sample_names) / columns)

    with PdfPages(output_file) as pdf:
        for index, mirna_id in enumerate(mirna_ids, start=1):
            filename = mirna_id.replace("*", "_star") + ".txt"
            figure, axes = plt.subplots(
                rows, columns, squeeze=False,
                figsize=(5.0 * columns, 4.8 * rows),
            )
            for sample_index, sample in enumerate(sample_names):
                axis = axes.flat[sample_index]
                total, matrix = _read_matrix(profile_root / sample / filename)
                _plot_sample(
                    axis, matrix, sample, total,
                    COLORS[sample_index % len(COLORS)],
                )
            for unused in axes.flat[len(sample_names):]:
                unused.set_visible(False)
            figure.suptitle(mirna_id, fontsize=13, fontweight="bold")
            figure.tight_layout(rect=(0, 0, 1, 0.96))
            pdf.savefig(figure)
            plt.close(figure)
            if index == 1 or index % 50 == 0 or index == len(mirna_ids):
                logging.info("Bubble plot pages: %d/%d", index, len(mirna_ids))
    logging.info("Bubble plot complete: %s", output_file)
