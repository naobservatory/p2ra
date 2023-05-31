#!/usr/bin/env python3
from pathlib import Path

import stats
from mgs import BioProject, Enrichment, MGSData
from pathogen_properties import NAType
from pathogens import pathogens

bioprojects = {
    "crits_christoph": BioProject("PRJNA661613"),
    "rothman": BioProject("PRJNA729801"),
    "spurbeck": BioProject("PRJNA924011"),
}


def start():
    fdir = Path("fits")
    fdir.mkdir(exist_ok=True)
    figdir = Path("fig")
    figdir.mkdir(exist_ok=True)
    mgs_data = MGSData.from_repo()
    predictor = "incidence"
    with open("estimate_table.tsv", "w") as f:
        print(
            "pathogen_name",
            "study",
            "fine_location",
            "state",
            "county",
            "date",
            "estimated",
            sep="\t",
            file=f,
        )
        for pathogen_name, pathogen in pathogens.items():
            if pathogen.pathogen_chars.na_type != NAType.RNA:
                continue
            for study, bioproject in bioprojects.items():
                samples = mgs_data.sample_attributes(
                    bioproject, enrichment=Enrichment.VIRAL
                )
                predictors: list[stats.Predictor]
                if predictor == "incidence":
                    predictors = pathogen.estimate_incidences()
                elif predictor == "prevalence":
                    predictors = pathogen.estimate_prevalences()
                else:
                    raise ValueError(
                        f"{predictor} must be one of 'incidence' or 'prevalence'"
                    )
                for attrs in samples.values():
                    try:
                        _ = stats.lookup_variable(attrs, predictors)
                    except AssertionError:
                        found = False
                    else:
                        found = True
                    print(
                        pathogen_name,
                        study,
                        attrs.fine_location,
                        attrs.state,
                        attrs.county,
                        attrs.date,
                        found,
                        sep="\t",
                        file=f,
                    )


if __name__ == "__main__":
    start()
