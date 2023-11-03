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
    TARGET_INCIDENCE = 0.01
    TARGET_THRESHOLD = 100
    viruses = [
        ("Norovirus (GII)", "incidence"),
        ("SARS-COV-2", "incidence"),
    ]
    study_labels = {
        "rothman": "Rothman",
        # "crits_christoph": "Crits-Christoph",
        "spurbeck": "Spurbeck",
    }

    for virus, predictor_type in viruses:
        for study in study_labels.keys():
            stats = data[virus, predictor_type, study, "Overall"]
            median = stats.percentiles[50]
            lower = stats.percentiles[25]
            upper = stats.percentiles[75]
            cumulative_incidence = np.logspace(-4, -1, 100)

            for value, value_str in [
                (median, "median"),
                (lower, "lower"),
                (upper, "upper"),
            ]:
                detection = round(
                    TARGET_THRESHOLD / (100 * value * TARGET_INCIDENCE)
                )
                print(
                    f"{virus} in {study} requires {value_str} {detection} reads for detection with a threshhold of {TARGET_THRESHOLD}  at {TARGET_INCIDENCE} incidence"
                )


if __name__ == "__main__":
    start()
