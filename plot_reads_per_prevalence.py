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
        ("HIV", "prevalence"),
    ]
    fig, axes = plt.subplots(1, len(viruses), sharey=True, figsize=(9, 2.75))
    for (virus, predictor_type), ax in zip(viruses, axes):
        stats = data[virus, predictor_type, "rothman", "Overall"]
        median = stats.percentiles[50]
        lower = stats.percentiles[25]
        upper = stats.percentiles[75]
        incidence = np.logspace(-4, -1, 100)
        for i, detection_threshold in enumerate([200, 20, 2]):
            # 0.01 = conversion from per 1% to true incidence
            scale_factor = detection_threshold * 0.01 / incidence
            color = f"C{i}"
            ax.loglog(
                incidence,
                scale_factor / median,
                color=color,
                label=f"{detection_threshold} reads",
            )
            ax.fill_between(
                incidence,
                scale_factor / upper,
                scale_factor / lower,
                color=color,
                alpha=0.2,
            )
        ax.set_title(virus)
        ax.set_xlabel("Prevalence")
        ax.grid()
        ax.set_xticks([1e-4, 1e-3, 1e-2, 1e-1])
    axes[0].set_ylabel("Reads required for detection")
    axes[-1].legend(
        title="Detection threshold", frameon=True, facecolor="w", framealpha=1
    )
    fig.savefig("reads_per_prevalence.png", bbox_inches="tight")


if __name__ == "__main__":
    start()
