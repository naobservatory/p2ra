#!/usr/bin/env python3

import csv
from dataclasses import dataclass

import matplotlib.pyplot as plt  # type: ignore

import numpy as np

PERCENTILES = [5, 25, 50, 75, 95]


@dataclass
class SummaryStats:
    mean: float
    std: float
    min: float
    percentiles: dict[int, float]
    max: float


def read_data() -> dict[tuple[str, str, str, str], SummaryStats]:
    data = {}
    with open("fits_summary.tsv") as datafile:
        reader = csv.DictReader(datafile, delimiter="\t")
        for row in reader:
            virus = row["tidy_name"]
            predictor_type = row["predictor_type"]
            study = row["study"]
            location = row["location"]
            data[virus, predictor_type, study, location] = SummaryStats(
                mean=float(row["mean"]),
                std=float(row["std"]),
                min=float(row["min"]),
                percentiles={p: float(row[f"{p}%"]) for p in PERCENTILES},
                max=float(row["max"]),
            )
    return data


def start():
    data = read_data()
    viruses = [
        ("Norovirus (GII)", "incidence"),
        ("SARS-COV-2", "incidence"),
        # ("HIV", "prevalence"),
    ]
    study_labels = {
        "rothman": "Rothman",
        "crits_christoph": "Crits-Christoph",
        "spurbeck": "Spurbeck",
    }
    n_viruses = len(viruses)
    # n_studies = len(study_labels)
    fig, (top_axes, bottom_axes) = plt.subplots(
        2,
        len(viruses),
        sharey=True,
        figsize=(3 * n_viruses, 5.5),
    )

    for axes in top_axes, bottom_axes:
        for (virus, predictor_type), ax in zip(viruses, axes):
            for i, study in enumerate(study_labels.keys()):
                stats = data[virus, predictor_type, study, "Overall"]
                median = stats.percentiles[50]
                lower = stats.percentiles[25]
                upper = stats.percentiles[75]
                cumulative_incidence = np.logspace(-4, -1, 100)
                detection_threshold = 1000 if axes is top_axes else 100
                fig.text(
                    0.5,
                    0.94 if axes is top_axes else 0.47,
                    f"Detection Threshold: {detection_threshold}",
                    ha="center",
                    va="center",
                    fontsize=10,
                )

                color = f"C{i}"
                ax.set_box_aspect(0.9)

                ax.loglog(
                    cumulative_incidence,
                    detection_threshold
                    / (100 * median * cumulative_incidence),
                    color=color,
                    label=f"{study_labels[study]}",
                )

                ax.fill_between(
                    cumulative_incidence,
                    detection_threshold / (100 * lower * cumulative_incidence),
                    detection_threshold / (100 * upper * cumulative_incidence),
                    color=color,
                    alpha=0.2,
                )
                if axes is top_axes:
                    title_text = ax.text(
                        0.5,
                        1.05,
                        virus,
                        transform=ax.transAxes,
                        horizontalalignment="center",
                    )

                ax.grid()
                ax.set_xticks([1e-4, 1e-3, 1e-2, 1e-1])
                ax.set_xticklabels(["0.01%", "0.1%", "1%", "10%"])
                ax.set_yticks([1e3, 1e6, 1e9, 1e12, 1e15])
                ax.set_xlim(1e-4, 1e-1)

    fig.subplots_adjust(hspace=0.4, wspace=0.2)

    for i, (top_ax, bottom_ax) in enumerate(zip(top_axes, bottom_axes)):
        if i == 0:
            for ax in top_ax, bottom_ax:
                ax.set_ylabel("Reads required for detection")
        else:
            for ax in top_ax, bottom_ax:
                ax.tick_params(axis="y", which="both", left=False, right=False)

        bottom_ax.set_xlabel("Cumulative Incidence")

    top_ax_title = fig.axes[0].set_title(
        "a",
        fontweight="bold",
    )
    bottom_ax_title = fig.axes[2].set_title(
        "b",
        fontweight="bold",
    )
    for ax in top_ax_title, bottom_ax_title:
        ax.set_position((-0.35, 1.05))
    legend = axes[1].legend(
        loc=(-1.2, -0.42),
        frameon=True,
        facecolor="w",
        framealpha=1,
        ncol=3,
    )

    fig.tight_layout
    fig.savefig("reads_per_incidence.png", bbox_inches="tight", dpi=600)


if __name__ == "__main__":
    start()
