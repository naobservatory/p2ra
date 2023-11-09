#!/usr/bin/env python3

import csv
from dataclasses import dataclass
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


def get_cost(
    reads_required=list, detection_threshold=int, TARGET_INCIDENCE=float
):
    DOLLAR_PER_1B_READS = 8000
    weeks_per_year = 52

    costs = []

    for reads in reads_required:
        cost = round(
            detection_threshold
            / (100 * reads * TARGET_INCIDENCE)
            * weeks_per_year
            * DOLLAR_PER_1B_READS
            / 1e9
        )

        tidy_cost = f"${cost:,}"
        costs.append(tidy_cost)
    tidy_median_cost, tidy_lower_cost, tidy_upper_cost = costs

    return tidy_median_cost, tidy_lower_cost, tidy_upper_cost


def start():
    data = read_data()
    TARGET_INCIDENCE = 0.01
    TARGET_THRESHOLDS = [10, 100, 1000]

    viruses = ["Norovirus (GII)", "SARS-COV-2"]
    study_labels = {
        "crits_christoph": "Crits-Christoph",
        "rothman": "Rothman",
        "spurbeck": "Spurbeck",
    }
    with open("cost_estimates.tsv", mode="w", newline="") as file:
        tsv_writer = csv.writer(file, delimiter="\t")
        tsv_writer.writerow(
            [
                "Virus",
                "Study",
                "50th %",
                "25th %",
                "75th %",
                "Detection Threshold",
            ]
        )
        for detection_threshold in TARGET_THRESHOLDS:
            for virus in viruses:
                geomean_dict = {
                    "median": [],
                    "lower": [],
                    "upper": [],
                }
                studies = study_labels.keys()

                for i, study in enumerate(studies):
                    stats = data[virus, "incidence", study, "Overall"]
                    study_median = stats.percentiles[50]
                    study_lower = stats.percentiles[25]
                    study_upper = stats.percentiles[75]

                    geomean_dict["median"].append(study_median)
                    geomean_dict["lower"].append(study_lower)
                    geomean_dict["upper"].append(study_upper)

                    (
                        tidy_median_cost,
                        tidy_lower_cost,
                        tidy_upper_cost,
                    ) = get_cost(
                        [study_median, study_lower, study_upper],
                        detection_threshold,
                        TARGET_INCIDENCE,
                    )

                    tsv_writer.writerow(
                        [
                            virus,
                            study_labels[study],
                            tidy_median_cost,
                            tidy_lower_cost,
                            tidy_upper_cost,
                            detection_threshold,
                        ]
                    )
                    if i == len(studies) - 1:
                        geomean_median = gmean(geomean_dict["median"])
                        geomean_lower = gmean(geomean_dict["lower"])
                        geomean_upper = gmean(geomean_dict["upper"])

                        (
                            tidy_median_cost,
                            tidy_lower_cost,
                            tidy_upper_cost,
                        ) = get_cost(
                            [geomean_median, geomean_lower, geomean_upper],
                            detection_threshold,
                            TARGET_INCIDENCE,
                        )

                        tsv_writer.writerow(
                            [
                                virus,
                                "Mean (geometric)",
                                tidy_median_cost,
                                tidy_lower_cost,
                                tidy_upper_cost,
                                detection_threshold,
                            ]
                        )


if __name__ == "__main__":
    start()
