#!/usr/bin/env python3

import csv
from dataclasses import dataclass

import matplotlib.pyplot as plt  # type: ignore
from matplotlib.lines import Line2D  # type: ignore
from scipy.stats import gmean
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


def plot_lines(
    ax: plt.Axes,
    median: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    label: str,
    color: str,
    linestyle: str,
    cumulative_incidence: int,
) -> None:
    ax.loglog(
        cumulative_incidence,
        median,
        color=color,
        label=label,
        linestyle=linestyle,
    )

    ax.fill_between(
        cumulative_incidence,
        lower,
        upper,
        color=color,
        alpha=0.2,
    )


def get_reads_required(
    data=dict,
    cumulative_incidence=int,
    detection_threshold=np.ndarray,
    virus=str,
    predictor_type=str,
    study=str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    stats = data[virus, predictor_type, study, "Overall"]

    median_reads = detection_threshold / (
        100 * stats.percentiles[50] * cumulative_incidence
    )
    lower_reads = detection_threshold / (
        100 * stats.percentiles[25] * cumulative_incidence
    )
    upper_reads = detection_threshold / (
        100 * stats.percentiles[75] * cumulative_incidence
    )

    return median_reads, lower_reads, upper_reads


def start():
    data = read_data()

    viruses = ["Norovirus (GII)", "SARS-COV-2"]
    study_labels = {
        "crits_christoph": "Crits-Christoph",
        "rothman": "Rothman",
        "spurbeck": "Spurbeck",
    }
    DETECTION_THRESHOLDS = [10, 100, 1000]

    fig, (top_axes, bottom_axes) = plt.subplots(
        len(viruses),
        len(DETECTION_THRESHOLDS),
        sharey=True,
        figsize=(9, 6),
    )

    for axes, detection_threshold in zip(
        zip(top_axes, bottom_axes), DETECTION_THRESHOLDS
    ):
        for virus, ax in zip(viruses, axes):
            geomean_dict = {
                "median": [],
                "lower": [],
                "upper": [],
            }
            studies = study_labels.keys()
            for i, study in enumerate(studies):
                study_median, study_lower, study_upper = get_reads_required(
                    data,
                    cumulative_incidence=np.logspace(-4, -1, 100),
                    detection_threshold=detection_threshold,
                    virus=virus,
                    predictor_type="incidence",
                    study=study,
                )

                geomean_dict["median"].append(study_median)
                geomean_dict["lower"].append(study_lower)
                geomean_dict["upper"].append(study_upper)

                cumulative_incidence = np.logspace(-4, -1, 100)
                detection_threshold = detection_threshold

                if virus == "Norovirus (GII)":
                    ax.set_title(
                        f"Detection Threshold: {detection_threshold}",
                        loc="center",
                    )

                color = f"C{i}"

                plot_lines(
                    ax=ax,
                    median=study_median,
                    lower=study_lower,
                    upper=study_upper,
                    label=f"{study_labels[study]}",
                    linestyle="-",
                    color=color,
                    cumulative_incidence=cumulative_incidence,
                )

                if i == len(studies) - 1:
                    geomean_median = gmean(geomean_dict["median"])
                    geomean_lower = gmean(geomean_dict["lower"])
                    geomean_upper = gmean(geomean_dict["upper"])
                    color = f"C{i + 1}"

                    plot_lines(
                        ax,
                        median=geomean_median,
                        lower=geomean_lower,
                        upper=geomean_upper,
                        label="Mean (geometric)",
                        linestyle="-",
                        color=color,
                        cumulative_incidence=cumulative_incidence,
                    )

                ax.set_xticks([1e-4, 1e-3, 1e-2, 1e-1])
                ax.set_xticklabels(["0.01%", "0.1%", "1%", "10%"])
                ax.set_yticks([1e3, 1e6, 1e9, 1e12, 1e15])
                ax.set_xlim(1e-4, 1e-1)
                ax.grid(
                    which="major",
                    linestyle="-",
                    linewidth=0.5,
                    color="gray",
                    alpha=0.7,
                )

    fig.subplots_adjust(hspace=0.4, wspace=0.2)

    for i, (top_ax, bottom_ax) in enumerate(zip(top_axes, bottom_axes)):
        if i == 0:
            for ax in top_ax, bottom_ax:
                ax.set_ylabel("Reads required for detection")
        else:
            for ax in top_ax, bottom_ax:
                ax.tick_params(axis="y", which="both", left=False, right=False)

        bottom_ax.set_xlabel("Cumulative Incidence")

    fig.axes[0].text(
        -0.35,
        1.08,
        "a",
        fontweight="bold",
        fontdict={"fontsize": 12},
        transform=fig.axes[0].transAxes,
    )

    fig.axes[3].text(
        -0.35,
        1.08,
        "b",
        fontweight="bold",
        fontdict={"fontsize": 12},
        transform=fig.axes[3].transAxes,
    )

    legend = fig.axes[4].legend(
        bbox_to_anchor=(0.5, -0.45),
        loc="lower center",
        ncol=4,
    )

    for ax in fig.axes:
        ax.tick_params(axis="x", which="minor", bottom=False)

    fig.tight_layout
    fig.show()
    fig.savefig("geomean_reads_required.png", bbox_inches="tight", dpi=600)


if __name__ == "__main__":
    start()
