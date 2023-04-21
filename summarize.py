#!/usr/bin/env python3

import pathogens


def start():
    for pathogen_name, pathogen in pathogens.pathogens.items():
        print(pathogen_name)
        for estimate in pathogen.estimate_prevalences():
            print("  %.2f per 100k" % estimate.infections_per_100k)


if __name__ == "__main__":
    start()
