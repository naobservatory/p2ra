import math

import csv
def unweighted_geom_mean_pos_rate_total_cases(
    pos_rate: float, total_cases: int
) -> float:
    return (pos_rate**0.7 * total_cases**0.3)


with open("combined_USA.tsv") as inf:
    head = None
    for row in csv.reader(inf, delimiter='\t'):
        if head is None:
            head = row
            continue
        percent_positive = float(row[head.index("PCR Detection")]) / 100
        total_cases = int(row[head.index("PCR_Detections")])
        total_tests = (100 / percent_positive) * total_cases
        date = row[head.index('RepWeekDate')]



        prevalence = float(unweighted_geom_mean_pos_rate_total_cases
        (percent_positive, total_cases))
        # Print all values with tabs in between

        print(
            date,
            prevalence,
            total_cases,
            percent_positive,
            total_tests,
            sep='\t'
        )
