#!/usr/bin/env python3

import csv
from dataclasses import dataclass
from scipy.stats import gmean  # type: ignore

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


def tidy_number(reads_required=int) -> str:
    sci_notation = f"{reads_required:.2e}"

    coefficient, exponent = sci_notation.split("e")

    # Remove the leading '+' from the exponent
    exponent = exponent.replace("+", "")
    # Remove the leading zero from the exponent if it's there
    if exponent.startswith("0") and len(exponent) > 1:
        exponent = exponent[1:]

    # Now replace the digits with superscript characters
    exponent = (
        exponent.replace("0", "⁰")
        .replace("1", "¹")
        .replace("2", "²")
        .replace("3", "³")
        .replace("4", "⁴")
        .replace("5", "⁵")
        .replace("6", "⁶")
        .replace("7", "⁷")
        .replace("8", "⁸")
        .replace("9", "⁹")
    )

    return f"{coefficient} x 10{exponent}"


# def tidy_number(reads_required=int) -> str:
#     scientific_notation = "{:.2e} reads".format(round(reads_required))


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
    with open("read_estimates.tsv", mode="w", newline="") as file:
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
                    (
                        study_median,
                        study_lower,
                        study_upper,
                    ) = get_reads_required(
                        data,
                        cumulative_incidence=TARGET_INCIDENCE,
                        detection_threshold=detection_threshold,
                        virus=virus,
                        predictor_type="incidence",
                        study=study,
                    )

                    geomean_dict["median"].append(study_median)
                    geomean_dict["lower"].append(study_lower)
                    geomean_dict["upper"].append(study_upper)

                    tidy_median = tidy_number(study_median)
                    tidy_lower = tidy_number(study_lower)
                    tidy_upper = tidy_number(study_upper)
                    tsv_writer.writerow(
                        [
                            virus,
                            study_labels[study],
                            tidy_median,
                            tidy_lower,
                            tidy_upper,
                            detection_threshold,
                        ]
                    )
                    if i == len(studies) - 1:
                        geomean_median = gmean(geomean_dict["median"])
                        geomean_lower = gmean(geomean_dict["lower"])
                        geomean_upper = gmean(geomean_dict["upper"])

                        tidy_median = tidy_number(geomean_median)
                        tidy_lower = tidy_number(geomean_lower)
                        tidy_upper = tidy_number(geomean_upper)

                        tsv_writer.writerow(
                            [
                                virus,
                                "Mean (geometric)",
                                tidy_median,
                                tidy_lower,
                                tidy_upper,
                                detection_threshold,
                            ]
                        )


if __name__ == "__main__":
    start()
