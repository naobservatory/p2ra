#!/usr/bin/env python3

import pathogens


def start():
    for pathogen_name, pathogen in pathogens.pathogens.items():
        print(pathogen_name)
        for estimate_name, estimate in pathogen.estimate_prevalences().items():
            print(
                "  %s: %.2f per 100k"
                % (estimate_name, estimate.infections_per_100k)
            )


if __name__ == "__main__":
    start()
