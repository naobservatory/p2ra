#!/usr/bin/env python3
import calendar
import sys

import pathogens
from pathogen_properties import Variable

MAX_ESTIMATES_FOR_PATHOGEN = 5


def pretty_date(estimate: Variable) -> str:
    start_date, end_date = estimate.get_dates()
    _, end_date_month_last_day = calendar.monthrange(
        end_date.year, end_date.month
    )
    if start_date == end_date:
        return str(start_date)
    elif start_date.year != end_date.year:
        return f"{start_date.year} to {end_date.year}"
    elif (
        start_date.month == 1
        and start_date.day == 1
        and end_date.month == 12
        and end_date.day == 31
    ):
        return str(start_date.year)
    elif (
        start_date.month == end_date.month
        and start_date.day == 1
        and end_date.day == end_date_month_last_day
    ):
        return f"{start_date.year}-{start_date.month:02d}"
    else:
        return f"{start_date} to {end_date}"


def start(pathogen_names):
    maximum = None if pathogen_names else MAX_ESTIMATES_FOR_PATHOGEN
    for pathogen_name, pathogen in pathogens.pathogens.items():
        if pathogen_names and pathogen_name not in pathogen_names:
            continue

        print(pathogen_name)

        to_print = []  # key, line

        def save(estimate: Variable, details: str):
            taxid = ""
            if estimate.taxid:
                taxid = "; %s" % estimate.taxid
            date = pretty_date(estimate)
            location = estimate.summarize_location()
            line = "%s (%s; %s%s)" % (
                details.rjust(18),
                location,
                date,
                taxid,
            )
            to_print.append(((location, date), line))

        for estimate in pathogen.estimate_prevalences():
            save(estimate, "%.2f per 100k" % estimate.infections_per_100k)

        for estimate in pathogen.estimate_incidences():
            save(
                estimate,
                "%.2f per 100k/y" % estimate.annual_infections_per_100k,
            )

        to_print.sort()
        for _, line in to_print[:maximum]:
            print(line)
        skipped = 0 if not maximum else len(to_print[maximum:])
        if skipped:
            print(("+ %s more" % skipped).rjust(18))


if __name__ == "__main__":
    start(sys.argv[1:])
