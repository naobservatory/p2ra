#!/usr/bin/env python3

import csv
from dataclasses import dataclass


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
    TARGET_INCIDENCE = [0.01, 0.001
    TARGET_THRESHOLDS = [100, 1000]

    DOLLAR_PER_1B_READS = 8000
    weeks_per_year = 52
    viruses = [
        ("Norovirus (GII)", "incidence"),
        ("SARS-COV-2", "incidence"),
    ]
    study_labels = {
        "rothman": "Rothman",
        "crits_christoph": "Crits-Christoph",
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
        for threshold in TARGET_THRESHOLDS:
            for virus, predictor_type in viruses:
                for study in study_labels.keys():
                    stats = data[virus, predictor_type, study, "Overall"]
                    median = stats.percentiles[50]
                    lower = stats.percentiles[25]
                    upper = stats.percentiles[75]
                    cumulative_incidence = np.logspace(-4, -1, 100)
                    median_cost = round(
                        threshold
                        / (100 * median * TARGET_INCIDENCE)
                        * weeks_per_year
                        * DOLLAR_PER_1B_READS
                        / 1e9
                    )
                    lower_cost = round(
                        threshold
                        / (100 * lower * TARGET_INCIDENCE)
                        * weeks_per_year
                        * DOLLAR_PER_1B_READS
                        / 1e9
                    )
                    upper_cost = round(
                        threshold
                        / (100 * upper * TARGET_INCIDENCE)
                        * weeks_per_year
                        * DOLLAR_PER_1B_READS
                        / 1e9
                    )

                    tidy_median_cost = f"${median_cost} (50%)"
                    tidy_lower_cost = f"${lower_cost} (25%)"
                    tidy_upper_cost = f"${upper_cost} (75%)"

                    # Write each estimate as a new row in the TSV.
                    tsv_writer.writerow(
                        [
                            virus,
                            study_labels[study],
                            tidy_median_cost,
                            tidy_lower_cost,
                            tidy_upper_cost,
                            threshold,
                        ]
                    )


if __name__ == "__main__":
    start()
