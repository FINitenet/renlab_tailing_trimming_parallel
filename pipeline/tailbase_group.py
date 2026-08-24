import logging
import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BASE_COLORS = {
    "A": "#00FFFF",
    "C": "#FF0000",
    "G": "#00FF00",
    "T": "#0000FF",
}
STACK_ORDER = ("T", "G", "C", "A")
LEGEND_ORDER = ("Tailing A", "Tailing C", "Tailing G", "Tailing T", "Trimming")


def _safe_name(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_")


def _read_metadata(path):
    metadata = pd.read_csv(path, sep="\t")
    if "sample_raw" not in metadata.columns and "sample" in metadata.columns:
        metadata = metadata.rename(columns={"sample": "sample_raw"})
    required = {"sample_raw", "tissue", "genotype", "replicate"}
    missing = sorted(required - set(metadata.columns))
    if missing:
        raise ValueError(
            "Group metadata is missing required columns: " + ", ".join(missing)
        )
    if metadata["sample_raw"].duplicated().any():
        duplicated = metadata.loc[metadata["sample_raw"].duplicated(), "sample_raw"].tolist()
        raise ValueError(f"Duplicate samples in group metadata: {duplicated}")
    return metadata


def _collect(metadata, table_root, max_length):
    records = []
    missing_files = []
    for row in metadata.itertuples(index=False):
        sample = str(row.sample_raw)
        sample_dir = table_root / sample
        tail_file = sample_dir / f"{sample}_tail_nt_summary.xlsx"
        trim_file = sample_dir / f"{sample}_trim_summary.xlsx"
        for path in (tail_file, trim_file):
            if not path.is_file():
                missing_files.append(str(path))
        if not tail_file.is_file() or not trim_file.is_file():
            continue

        tail = pd.read_excel(tail_file)
        trim = pd.read_excel(trim_file)
        tail_required = {"tail_N", "A_percentage", "C_percentage", "G_percentage", "T_percentage"}
        trim_required = {"trim_N", "percentage"}
        if not tail_required.issubset(tail.columns):
            missing = sorted(tail_required - set(tail.columns))
            raise ValueError(f"{tail_file} is missing columns: {missing}")
        if not trim_required.issubset(trim.columns):
            missing = sorted(trim_required - set(trim.columns))
            raise ValueError(f"{trim_file} is missing columns: {missing}")

        common = {
            "sample_raw": sample,
            "tissue": row.tissue,
            "genotype": row.genotype,
            "replicate": row.replicate,
        }
        for tail_row in tail[tail["tail_N"].between(1, max_length)].itertuples(index=False):
            for base in ("A", "C", "G", "T"):
                records.append({
                    **common,
                    "event": "Tailing",
                    "base": base,
                    "length": int(tail_row.tail_N),
                    "percentage": float(getattr(tail_row, f"{base}_percentage")),
                })
        for trim_row in trim[trim["trim_N"].between(1, max_length)].itertuples(index=False):
            records.append({
                **common,
                "event": "Trimming",
                "base": "Trimming",
                "length": int(trim_row.trim_N),
                "percentage": float(trim_row.percentage),
            })

    if missing_files:
        preview = "\n".join(missing_files[:10])
        suffix = "" if len(missing_files) <= 10 else f"\n... and {len(missing_files) - 10} more"
        raise FileNotFoundError(f"Missing tailbase tables:\n{preview}{suffix}")
    return pd.DataFrame.from_records(records)


def _plot_tissue(tissue, genotypes, summary, tail_total, max_length, output_dir):
    ncols = min(2, len(genotypes))
    nrows = math.ceil(len(genotypes) / ncols)
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(5 * ncols, 3.5 * nrows),
        constrained_layout=True, sharex=True, sharey=True,
    )
    axes = np.atleast_1d(axes).ravel()
    x = np.arange(1, max_length + 1)

    for axis, genotype in zip(axes, genotypes):
        subset = summary[(summary["tissue"] == tissue) & (summary["genotype"] == genotype)]
        bottom = np.zeros(len(x))
        for base in STACK_ORDER:
            values = subset[(subset["event"] == "Tailing") & (subset["base"] == base)]
            means = values.set_index("length")["mean"].reindex(x, fill_value=0).to_numpy()
            axis.bar(
                x, means, bottom=bottom, width=0.78,
                color=BASE_COLORS[base], alpha=0.9, label=f"Tailing {base}",
            )
            bottom += means

        total = tail_total[
            (tail_total["tissue"] == tissue) & (tail_total["genotype"] == genotype)
        ].set_index("length")
        axis.errorbar(
            x,
            total["mean"].reindex(x, fill_value=0).to_numpy(),
            yerr=total["sd"].reindex(x, fill_value=0).fillna(0).to_numpy(),
            fmt="none", ecolor="black", elinewidth=0.8, capsize=2, zorder=5,
        )

        trimming = subset[subset["event"] == "Trimming"].set_index("length")
        axis.bar(
            x,
            -trimming["mean"].reindex(x, fill_value=0).to_numpy(),
            yerr=trimming["sd"].reindex(x, fill_value=0).fillna(0).to_numpy(),
            width=0.78, color="#FFA500", alpha=0.9, label="Trimming", capsize=2,
        )
        axis.set_title(str(genotype), fontsize=10)
        axis.set_xlabel("Modification length (nt)")
        axis.set_ylabel("Percentage (%)\n(trimming shown below zero)")
        axis.spines[["top", "right"]].set_visible(False)
        axis.tick_params(labelsize=8)

    for axis in axes[len(genotypes):]:
        axis.set_visible(False)

    handles, labels = axes[0].get_legend_handles_labels()
    ordered = [(handles[labels.index(label)], label) for label in LEGEND_ORDER]
    axes[min(1, len(genotypes) - 1)].legend(
        [item[0] for item in ordered], [item[1] for item in ordered],
        frameon=False, fontsize=8,
    )
    fig.suptitle(f"{tissue}: miRNA tailing and trimming", fontsize=13)
    stem = output_dir / f"{_safe_name(tissue)}_miRNA_tailing_trimming"
    fig.savefig(stem.with_suffix(".pdf"), metadata={"CreationDate": None})
    fig.savefig(stem.with_suffix(".png"), dpi=300)
    plt.close(fig)


def run(config):
    metadata_path = Path(config["tailbase_group_metadata"])
    table_root = Path(config["tailbase_doc_dir"])
    output_dir = Path(config["tailbase_group_output"])
    max_length = int(config["tailbase_group_max_length"])

    if config["dry_run"]:
        logging.info(
            "Would create nucleotide-resolved tailbase group summaries: %s -> %s",
            metadata_path, output_dir,
        )
        return

    metadata = _read_metadata(metadata_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    long_table = _collect(metadata, table_root, max_length)
    long_table.to_csv(output_dir / "miRNA_tailing_trimming_long.tsv", sep="\t", index=False)

    summary = (
        long_table.groupby(["tissue", "genotype", "event", "base", "length"], sort=False)["percentage"]
        .agg(mean="mean", sd="std")
        .reset_index()
    )
    summary.to_csv(output_dir / "miRNA_tailing_trimming_group_summary.tsv", sep="\t", index=False)

    tail_total = (
        long_table[long_table["event"] == "Tailing"]
        .groupby(["sample_raw", "tissue", "genotype", "replicate", "length"], as_index=False)["percentage"]
        .sum()
        .groupby(["tissue", "genotype", "length"], sort=False)["percentage"]
        .agg(mean="mean", sd="std")
        .reset_index()
    )
    tail_total.to_csv(output_dir / "miRNA_tailing_total_group_summary.tsv", sep="\t", index=False)

    tissues = list(dict.fromkeys(metadata["tissue"].tolist()))
    for tissue in tissues:
        tissue_rows = metadata[metadata["tissue"] == tissue]
        genotypes = list(dict.fromkeys(tissue_rows["genotype"].tolist()))
        _plot_tissue(tissue, genotypes, summary, tail_total, max_length, output_dir)
    logging.info(
        "Tailbase group summary contains %d samples, %d tissues, and %d rows",
        metadata["sample_raw"].nunique(), len(tissues), len(long_table),
    )
