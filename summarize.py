#!/usr/bin/env python3
import sys

import pathogens


def start(pathogen_names):
    for pathogen_name, pathogen in pathogens.pathogens.items():
        if pathogen_names and pathogen_name not in pathogen_names:
            continue

        print(pathogen_name)
        for estimate in pathogen.estimate_prevalences():
            print(
                "  %.2f per 100k (%s; %s)"
                % (
                    estimate.infections_per_100k,
                    estimate.summarize_location(),
                    estimate.summarize_date(),
                )
            )


if __name__ == "__main__":
    start(sys.argv[1:])
