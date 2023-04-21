#!/usr/bin/env python3
import sys

import pathogens

MAX_ESTIMATES_FOR_PATHOGEN = 5


def start(pathogen_names):
    for pathogen_name, pathogen in pathogens.pathogens.items():
        if pathogen_names and pathogen_name not in pathogen_names:
            continue

        print(pathogen_name)
        for n, estimate in enumerate(pathogen.estimate_prevalences()):
            if n > MAX_ESTIMATES_FOR_PATHOGEN:
                print("  ...")
                break
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
