import csv

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
with open("2015-2016/ORXHPV_I.csv") as f:
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

for age in sorted(rates):
    positives = rates[age]["1"]
    negatives = rates[age]["2"]
    total = positives + negatives
    print("%s\t%s\t%s\t%s" % (age, positives, negatives, total))


