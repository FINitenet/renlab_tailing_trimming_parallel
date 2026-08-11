#!/usr/bin/env python3
"""Generate anonymous, deterministic figures used in the README.

These plots illustrate the shape of workflow outputs. They do not contain
experimental measurements or sample identifiers from a real project.
"""

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/renlab-tailing-readme-matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


OUTPUT_DIR = Path(__file__).resolve().parent / "images"
RNA_TYPES = (
    "unassigned", "snoRNA", "snRNA", "rRNA", "tRNA", "lncRNA",
    "protein_coding", "transposable_element", "hc-siRNA",
    "phasi_tasiRNA", "miRNA_primary_transcript",
)
RNA_COLORS = {
    "unassigned": "#a6cee3", "snoRNA": "#33a02c", "snRNA": "#fb9a99",
    "rRNA": "#1f78b4", "tRNA": "#b2df8a", "lncRNA": "#e31a1c",
    "protein_coding": "#fdbf6f", "transposable_element": "#ff7f00",
    "hc-siRNA": "#cab2d6", "phasi_tasiRNA": "#6a3d9a",
    "miRNA_primary_transcript": "#ffff99",
}


def _finish(figure, filename):
    figure.text(
        0.995, 0.005, "Schematic example — not experimental data",
        ha="right", va="bottom", fontsize=7, color="#666666",
    )
    figure.savefig(OUTPUT_DIR / filename, dpi=180, bbox_inches="tight")
    plt.close(figure)


def qc_example():
    positions = np.arange(1, 31)
    quality = 35.8 - 0.035 * (positions - 10) ** 2
    lengths = np.arange(12, 31)
    abundance = np.exp(-0.5 * ((lengths - 21) / 2.2) ** 2)
    abundance = abundance / abundance.sum() * 100

    figure, axes = plt.subplots(1, 2, figsize=(10, 3.6))
    axes[0].plot(positions, quality, color="#4c78a8", linewidth=2)
    axes[0].axhspan(28, 40, color="#d8efcf", zorder=-1)
    axes[0].axhspan(20, 28, color="#fff0bd", zorder=-1)
    axes[0].set(xlabel="Position in read (nt)", ylabel="Phred quality", ylim=(20, 40))
    axes[0].set_title("Quality", fontweight="bold")
    axes[1].bar(lengths, abundance, color="#54a9a6", width=0.8)
    axes[1].set(xlabel="Read length (nt)", ylabel="Percentage (%)")
    axes[1].set_title("Length", fontweight="bold")
    for axis in axes:
        axis.grid(axis="y", color="#d0d0d0", linewidth=0.5)
        axis.set_axisbelow(True)
    figure.suptitle("Sample_A", fontweight="bold")
    figure.tight_layout(rect=(0, 0.03, 1, 0.94))
    _finish(figure, "qc_example.png")


def bubble_example():
    coordinates = np.arange(10, -1, -1)
    x_values, y_values = np.meshgrid(coordinates, coordinates)
    figure, axes = plt.subplots(1, 2, figsize=(8.6, 4))
    for index, (axis, sample, color) in enumerate(
        zip(axes, ("Sample_A", "Sample_B"), ("#9b4b8f", "#4c9a51"))
    ):
        matrix = np.zeros((11, 11), dtype=float)
        matrix[10, 10] = 72 - index * 8
        matrix[10, 9] = 12 + index * 3
        matrix[9, 10] = 7 + index * 2
        matrix[10, 8] = 4
        matrix[8, 10] = 3
        matrix[9, 9] = 2
        sizes = matrix.ravel() / matrix.sum() * 1350
        axis.scatter(x_values.ravel(), y_values.ravel(), s=sizes, color=color, edgecolors="none")
        axis.set(xlim=(11, -1), ylim=(-1, 11), xlabel="Length of trimming", ylabel="Length of tailing")
        axis.set_xticks(range(11))
        axis.set_yticks(range(11))
        axis.grid(color="#b8b8b8", linewidth=0.5)
        axis.set_aspect("equal")
        axis.set_title(f"{sample}\ntotal reads: {12480 + index * 1730}", fontweight="bold", fontsize=10)
    figure.suptitle("ath-miR-demo", fontweight="bold")
    figure.tight_layout(rect=(0, 0.04, 1, 0.94))
    _finish(figure, "bubble_example.png")


def length_example():
    lengths = np.arange(12, 31)
    whole = np.exp(-0.5 * ((lengths - 21.5) / 2.1) ** 2)
    gmc = np.exp(-0.5 * ((lengths - 20.5) / 1.8) ** 2)
    whole = whole / whole.sum() * 100
    gmc = gmc / gmc.sum() * 100
    figure, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    for axis, values, title, color in (
        (axes[0], whole, "Whole reads", "#4c78a8"),
        (axes[1], gmc, "5GMC", "#9b4b8f"),
    ):
        axis.bar(lengths, values, color=color, width=0.8, edgecolor="white", linewidth=0.4)
        axis.set_title(title, fontweight="bold")
        axis.set_xlabel("Length (nt)")
        axis.set_xticks(lengths[::2])
        axis.grid(axis="y", color="#d0d0d0", linewidth=0.5)
        axis.set_axisbelow(True)
    axes[0].set_ylabel("Percentage of mapped miRNA reads (%)")
    figure.suptitle("Sample_A", fontweight="bold")
    figure.tight_layout(rect=(0, 0.04, 1, 0.93))
    _finish(figure, "length_example.png")


def tailbase_example():
    lengths = np.arange(1, 7)
    bases = {
        "A": np.array([24.0, 11.0, 5.0, 2.6, 1.2, 0.7]),
        "T": np.array([11.0, 6.5, 3.2, 1.6, 0.8, 0.4]),
        "C": np.array([6.0, 3.2, 1.6, 0.8, 0.4, 0.2]),
        "G": np.array([4.0, 2.1, 1.1, 0.5, 0.25, 0.12]),
    }
    colors = {"A": "#00FFFF", "T": "#0000FF", "C": "#FF0000", "G": "#00FF00"}
    trimming = -np.array([32.0, 15.0, 7.0, 3.4, 1.6, 0.8])
    figure, axis = plt.subplots(figsize=(8.6, 4.4))
    bottom = np.zeros_like(lengths, dtype=float)
    for base in ("A", "T", "C", "G"):
        axis.bar(lengths, bases[base], bottom=bottom, label=f"Tailing: {base}", color=colors[base])
        bottom += bases[base]
    axis.bar(lengths, trimming, label="Trimming", color="#FFA500")
    axis.axhline(0, color="#333333", linewidth=0.7)
    axis.set(
        xlabel="Tailing and trimming length (nt)",
        ylabel="Relative Percentage (%)", xticks=lengths, ylim=(-100, 100),
    )
    axis.set_title("Sample_A", fontweight="bold")
    axis.legend(
        frameon=False, ncol=1, fontsize=8,
        bbox_to_anchor=(1.01, 1), loc="upper left",
    )
    axis.grid(axis="y", color="#d0d0d0", linewidth=0.5)
    axis.set_axisbelow(True)
    figure.tight_layout(rect=(0, 0.04, 1, 1))
    _finish(figure, "tailbase_example.png")


def rnatype_example():
    samples = ("Sample_A", "Sample_B", "Sample_C")
    overall = np.array([
        [7, 2, 1, 13, 6, 3, 5, 20, 18, 7, 18],
        [5, 1, 1, 10, 7, 4, 6, 17, 23, 8, 18],
        [8, 2, 1, 15, 5, 3, 5, 18, 16, 9, 18],
    ], dtype=float)
    lengths = np.arange(18, 29)
    by_length = np.empty((len(lengths), len(RNA_TYPES)))
    for row, length in enumerate(lengths):
        mirna = 7 + 38 * np.exp(-0.5 * ((length - 21) / 1.15) ** 2)
        hc = 9 + 26 * np.exp(-0.5 * ((length - 24) / 1.35) ** 2)
        other = np.array([7, 2, 1, 13, 7, 3, 5, 17, hc, 7, mirna], dtype=float)
        by_length[row] = other / other.sum() * 100

    figure, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
    for axis, labels, matrix, title in (
        (axes[0], samples, overall, None),
        (axes[1], lengths, by_length, "Sample_A"),
    ):
        bottom = np.zeros(matrix.shape[0])
        handles = []
        for category in reversed(RNA_TYPES):
            index = RNA_TYPES.index(category)
            bars = axis.bar(labels, matrix[:, index], bottom=bottom, color=RNA_COLORS[category], label=category)
            bottom += matrix[:, index]
            handles.append(bars)
        axis.set_ylim(0, 100)
        axis.set_ylabel("Percentage (%)")
        if title:
            axis.set_title(title, fontweight="bold")
        axis.grid(axis="y", color="#d0d0d0", linewidth=0.5)
        axis.set_axisbelow(True)
    axes[1].set_xlabel("Read length (nt)")
    handles, labels = axes[1].get_legend_handles_labels()
    axes[1].legend(handles[::-1], labels[::-1], title="RNA type", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False, fontsize=8)
    figure.tight_layout(rect=(0, 0.04, 1, 1))
    _finish(figure, "rnatype_example.png")


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    qc_example()
    bubble_example()
    length_example()
    tailbase_example()
    rnatype_example()


if __name__ == "__main__":
    main()
