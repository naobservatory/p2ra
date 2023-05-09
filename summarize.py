#!/usr/bin/env python3
import calendar
import sys

import pathogens

MAX_ESTIMATES_FOR_PATHOGEN = 5


def start(pathogen_names):
    maximum = None if pathogen_names else MAX_ESTIMATES_FOR_PATHOGEN
    for pathogen_name, pathogen in pathogens.pathogens.items():
        if pathogen_names and pathogen_name not in pathogen_names:
            continue

        print(pathogen_name)

        to_print = []  # key, line

        for n, estimate in enumerate(pathogen.estimate_prevalences()):
            date = "no date"
            date_summary = estimate.summarize_date()
            if date_summary:
                start_date, end_date = date_summary
                _, end_date_month_last_day = calendar.monthrange(
                    end_date.year, end_date.month
                )
                if start_date == end_date:
                    date = start_date
                elif (
                    start_date.month == 1
                    and start_date.day == 1
                    and end_date.month == 12
                    and end_date.day == 31
                ):
                    date = start_date.year
                elif (
                    start_date.month == end_date.month
                    and start_date.day == 1
                    and end_date.day == end_date_month_last_day
                ):
                    date = f"{start_date.year}-{start_date.month:02d}"
                else:
                    date = f"{start_date} to {end_date}"

            location = estimate.summarize_location()

            prevalence = "%.2f per 100k" % estimate.infections_per_100k
            taxid = ""
            if estimate.taxid:
                taxid = "; %s" % estimate.taxid
            line = "%s (%s; %s%s)" % (
                prevalence.rjust(18),
                location,
                date,
                taxid,
            )
            to_print.append(((location, str(date)), line))

        to_print.sort()
        for _, line in to_print[:maximum]:
            print(line)
        skipped = 0 if not maximum else len(to_print[maximum:])
        if skipped:
            print(("+ %s more" % skipped).rjust(18))


if __name__ == "__main__":
    start(sys.argv[1:])
