#!/usr/bin/env python3
import itertools
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
    mgs_data = MGSData.from_repo()

    missing = []

    for pathogen_name, pathogen in pathogens.items():
        if pathogen_name in ["west_nile_virus", "dengue"]:
            continue # We don't care about these
        
        if pathogen.pathogen_chars.na_type != NAType.RNA:
            continue
        
        incidences = pathogen.estimate_incidences()
        prevalences = pathogen.estimate_prevalences()

        for study, bioproject in bioprojects.items():
            for sample, sample_attributes in mgs_data.sample_attributes(
                bioproject, enrichment=Enrichment.VIRAL
            ).items():
                # Every sample should have at least one match.
                if not stats.lookup_variables(
                    sample_attributes, itertools.chain(incidences, prevalences)
                ):
                    missing.append(
                        (
                            study,
                            sample_attributes.date,
                            sample_attributes.fine_location,
                            pathogen_name,
                        )
                    )

    if not missing:
        return

    missing.sort()
    print("Missing:")
    for attrs in missing:
        print("  ", *attrs)
    exit(1)


if __name__ == "__main__":
    start()
