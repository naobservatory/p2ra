#!/usr/bin/env python3
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns  # type: ignore

sys.path.append("..")

from pathogens import pathogens

##########################################

# First, learn how to create Bayesian test statistics. Posteriors are in fits.tsv, but I don't know where the priors are.


#########################################
def nucleic_acid(pathogen: str) -> str:
    return pathogens[pathogen].pathogen_chars.na_type.value


def selection_round(pathogen: str) -> str:
    return pathogens[pathogen].pathogen_chars.selection.value


def study_name(study: str) -> str:
    return {
        "brinch": "Brinch (DNA)",
        "crits_christoph": "Crits-Christoph",
        "rothman": "Rothman",
        "spurbeck": "Spurbeck",
    }[study]


def filter_data(
    virus: str, study: str, location: str, data: pd.DataFrame
) -> pd.DataFrame:
    # Filter the data based on the given parameters
    filtered_data = data[
        (data["pathogen"] == virus)
        & (data["study"] == study_name(study))
        & (data["location"] == location)
    ]
    return filtered_data


def main():
    parent_dir = Path("..")
    input_df = pd.read_csv(parent_dir / "input.tsv", sep="\t")
    input_df["study"] = input_df.study.map(study_name)

    # Example: Filter data for a specific virus, study, and location
    virus = "SARS-COV-2"
    study = "spurbeck"
    location = "A"
    filtered_df = filter_data(virus, study, location, input_df)

    # Print the filtered data
    print(filtered_df)


if __name__ == "__main__":
    main()
