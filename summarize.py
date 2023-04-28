#!/usr/bin/env python3
import sys

import pathogens

MAX_ESTIMATES_FOR_PATHOGEN = 5


def start(pathogen_names):
    for pathogen_name, pathogen in pathogens.pathogens.items():
        if pathogen_names and pathogen_name not in pathogen_names:
            continue

        print(pathogen_name)
        skipped = 0
        for n, estimate in enumerate(pathogen.estimate_prevalences()):
            if n > MAX_ESTIMATES_FOR_PATHOGEN:
                skipped += 1
                continue

            date = "no date"
            date_summary = estimate.summarize_date()
            if date_summary:
                start_date, end_date = date_summary
                if start_date == end_date:
                    date = start_date
                elif start_date.year != end_date.year:
                    date = "%s to %s" % (start_date.year, end_date.year)
                elif (
                    start_date.month == 1
                    and start_date.day == 1
                    and end_date.month == 12
                    and end_date.day == 12
                ):
                    date = start_date.year
                elif start_date.month != end_date.month:
                    date = start_date.year
                else:
                    date = "%s to %s" % (start_date, end_date)

            print(
                "  %.2f per 100k (%s; %s)"
                % (
                    estimate.infections_per_100k,
                    estimate.summarize_location(),
                    date,
                )
            )
        if skipped:
            print("  + %s more" % skipped)


if __name__ == "__main__":
    start(sys.argv[1:])
