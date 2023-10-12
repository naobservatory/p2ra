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
        ("SARS-CoV-2", "incidence"),
        ("HIV", "prevalence"),
    ]
    study_labels = {
        "rothman": "Rothman",
        "crits_christoph": "Crits-Christoph",
        "spurbeck": "Spurbeck",
    }
    fig, axes = plt.subplots(1, len(viruses), sharey=True, figsize=(9, 2.75))
    for (virus, predictor_type), ax in zip(viruses, axes):
        for i, study in enumerate(study_labels.keys()):
            stats = data[virus, predictor_type, study, "Overall"]
            median = stats.percentiles[50]
            lower = stats.percentiles[25]
            upper = stats.percentiles[75]
            cumulative_incidence = np.logspace(-4, -1, 100)
            detection_threshold = 100

            color = f"C{i}"

            ax.loglog(
                cumulative_incidence,
                detection_threshold / (100 * median * cumulative_incidence),
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
            ax.set_title(virus)
            ax.set_xlabel("Cumulative Incidence")
            ax.grid()
            ax.set_xticks([1e-4, 1e-3, 1e-2, 1e-1])

    fig.subplots_adjust(top=0.8)
    for ax_number in 1, 2:
        axes[ax_number].tick_params(
            axis="y", which="both", left=False, right=False
        )

    plt.gca().invert_yaxis()

    axes[0].set_ylabel("Reads required for detection")
    legend = axes[1].legend(
        loc=(-0.5, -0.42),
        frameon=True,
        facecolor="w",
        framealpha=1,
        ncol=3,
    )

    fig.savefig("reads_per_prevalence.png", bbox_inches="tight", dpi=600)


if __name__ == "__main__":
    start()
