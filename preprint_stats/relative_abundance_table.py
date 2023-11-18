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


def tidy_number(reads_required=int) -> str:
    sci_notation = f"{reads_required:.2e}"

    coefficient, exponent = sci_notation.split("e")

    is_negative = exponent.startswith("-")
    if is_negative:
        exponent = exponent[1:]

    # Remove the leading zero from the exponent if it's there
    exponent = exponent.lstrip("0")

    # Add back the superscript minus if the exponent was negative
    if is_negative:
        exponent = "⁻" + exponent

    # Now replace the digits with superscript characters
    superscript_map = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
    exponent = exponent.translate(superscript_map)

    return f"{coefficient} x 10{exponent}"


def read_data() -> dict[tuple[str, str, str, str], SummaryStats]:
    data = {}
    with open("fits_summary.tsv") as datafile:
        reader = csv.DictReader(datafile, delimiter="\t")
        for row in reader:
            virus = row["tidy_name"]
            predictor_type = row["predictor_type"]
            study = row["study"]
            location = row["location"]
            if (
                virus == "AAV5"
            ):  # FIXME: Remove this when AAV5 is dropped earlier.
                continue
            data[virus, predictor_type, study, location] = SummaryStats(
                mean=tidy_number(float(row["mean"])),
                std=tidy_number(float(row["std"])),
                min=tidy_number(float(row["min"])),
                percentiles={
                    p: tidy_number(float(row[f"{p}%"])) for p in PERCENTILES
                },
                max=tidy_number(float(row["max"])),
            )
    return data


def create_tsv():
    data = read_data()
    viruses = set()
    for entry in data.keys():
        virus, predictor_type = entry[:2]
        viruses.add((virus, predictor_type))

    sorted_viruses = sorted(viruses, key=lambda x: (x[1], x[0]))

    study_tidy = {
        "rothman": "Rothman",
        "crits_christoph": "Crits-Christoph",
        "spurbeck": "Spurbeck",
        "brinch": "Brinch",
    }

    headers = ["Virus", "Study", "Median", "Lower", "Upper"]

    with open("output_summary.tsv", "w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=headers, delimiter="\t")
        writer.writeheader()

        for virus, predictor_type in sorted_viruses:
            studies = ["rothman", "crits_christoph", "spurbeck"] + (
                ["brinch"] if predictor_type == "prevalence" else []
            )
            for study in studies:
                stats = data[virus, predictor_type, study, "Overall"]
                writer.writerow(
                    {
                        "Virus": virus,
                        "Study": study_tidy[study],
                        "Median": stats.percentiles[50],
                        "Lower": stats.percentiles[5],
                        "Upper": stats.percentiles[95],
                    }
                )


def start():
    create_tsv()


if __name__ == "__main__":
    start()
