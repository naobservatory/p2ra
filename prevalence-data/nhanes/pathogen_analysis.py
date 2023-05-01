import csv

import os

from collections import Counter, defaultdict

demographics = {}

with open("2015-2016/DEMO_I.csv") as f:
    reader = csv.reader(f)
    cols = None
    for row in reader:
        if not cols:
            cols = row
        else:
            demographics[row[cols.index("seqn")]] = {
                "age": row[cols.index("ridageyr")],
            }

rates = defaultdict(Counter)

with open("2015-2016/ORHPV_I.csv") as f:
    reader = csv.reader(f)

    cols = None
    for row in reader:
        if not cols:
            cols = row
        else:
            overall = "3"
            for col, val in zip(cols, row):
                if col in ["seqn", "orxhpv", "orxgh", "orxgl"]:
                    continue
                if val == "1":
                    overall = val
                elif val == "2" and overall != "1":
                    overall = val
            if overall != row[cols.index("orxhpv")]:
                print(
                    "ERROR: %s %s %r"
                    % (overall, row[cols.index("orxhpv")], row)
                )
            rates[demographics[row[cols.index("seqn")]]["age"]][overall] += 1
# for age in sorted(rates):
#     positives = rates[age]["1"]
#     negatives = rates[age]["2"]
#     total = positives + negatives
#     print("%s\t%s\t%s\t%s" % (age, positives, negatives, total))


# rates = defaultdict(Counter)


HEPA_15_16_rates = defaultdict(Counter)


with open("2015-2016_variable_dict.tsv") as f:
    reader = csv.reader(f, delimiter="\t")
    cols = None
    # list of all variables in the 2015-2016 data, without "_I.csv"
    vars_2015 = [
        fname.removesuffix(".csv") for fname in os.listdir("2015-2016")
    ]
    print(vars_2015)
    with open("2015-2016_relevant_variables.tsv", "w") as f:
        writer = csv.writer(f, delimiter="\t")
        for row in reader:
            if not cols:
                cols = row
                writer.writerow(row)
            else:
                if row[cols.index("Data File Name")] in vars_2015:
                    writer.writerow(row)
        # Turn a list of lists into a csv file
    with open("2015-2016/HEPA_I.csv") as f:
        reader = csv.reader(f)
        cols = None
        for row in reader:
            if not cols:
                cols = row
                print(cols)
            else:
                hep_a_status = row[cols.index("lbxha")]

                HEPA_15_16_rates[demographics[row[cols.index("seqn")]]["age"]][
                    hep_a_status
                ] += 1


for age in sorted(HEPA_15_16_rates):
    positives = HEPA_15_16_rates[age]["1"]
    negatives = HEPA_15_16_rates[age]["2"]
    indeterminate = HEPA_15_16_rates[age]["3"]
    missing = HEPA_15_16_rates[age]["NA"]
    total = positives + negatives
    print("%s\t%s\t%s\t%s" % (age, positives, negatives, total)) # Could use zfill here to make sorting easier.

with open("2015-2016/HEPB_I.csv") as f:
    reader = csv.reader(f)
    cols = None
    for row in reader:
        if not cols:
            cols = row
        else:

