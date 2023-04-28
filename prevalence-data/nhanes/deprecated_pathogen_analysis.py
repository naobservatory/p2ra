import csv
import os
from collections import Counter, defaultdict


data_dictionary = {}
sum = 0

nhanes_cycles = ["2015-2016", "2017-2018"]

for cycle in nhanes_cycles:
    for f_name in os.listdir(f"{cycle}"):
        # Set the first part of the file name as a new key in the dictionary
        data_dictionary[f_name.split("_")[0]] = []
with open("2017-2018_variable_dict.tsv") as f:
    reader = csv.reader(f, delimiter="\t")
    keys = list(data_dictionary.keys())
    for row in reader:
        variable_abbrevation = row[2].split("_")[0]
        if variable_abbrevation in keys:
            data_dictionary[variable_abbrevation].append(row[1].strip())
            # remove trailing whitespace
            print("%s, $s" % (variable_abbrevation, row[1].strip()))











by the Aedes aegypti and Aedes albopictus mosquitoes, present in tropical and subtropical regions.








# TODO: Add second variable_dict for 2015 and 2016 for the variable description

# Pulling out rows in 2017-2018_variable_dict.tsv that match the variables

# with open("2017-2018_variable_dict.tsv") as f:
#     reader = csv.reader(f, delimiter="\t")
#     for key in data_dictionary:
#         for row in reader:
#             if code == row[0]:
#                 # Adding all column values in that row to the dictionary
#                 data_dictionary[code].append(row[1])

# data_dictionary = {code: [] for code in cols}

# with open("2017-2018_variable_dict.tsv") as f:
#     reader = csv.reader(f, delimiter="\t")
#     for row in reader:
#         code = row[0]
#         if code in data_dictionary:
#             # Adding all column values in that row to the dictionary
#             data_dictionary[code].append(row[1])


# print(data_dictionary)

# # demographics = {}

# # with open("data/sources/csv/2015-2016/DEMO_I.csv") as f:
# #     reader = csv.reader(f)
# #     cols = None
# #     for row in reader:
# #         if not cols:
# #             cols = row
# #         else:
# #             demographics[row[cols.index("seqn")]] = {
# #                 "age": row[cols.index("ridageyr")],
# #             }
# # rates = defaultdict(Counter)
# # with open("data/sources/csv/2015-2016/ORHPV_I.csv") as f:
# #     reader = csv.reader(f)

# #     cols = None
# #     for row in reader:
# #         if not cols:
# #             cols = row
# #         else:
# #             overall = "3"
# #             for col, val in zip(cols, row):
# #                 if col in ["seqn", "orxhpv", "orxgh", "orxgl"]:
# #                     continue
# #                 if val == "1":
# #                     overall = val
# #                 elif val == "2" and overall != "1":
# #                     overall = val
# #             if overall != row[cols.index("orxhpv")]:
# #                 print(
# #                     "ERROR: %s %s %r"
# #                     % (overall, row[cols.index("orxhpv")], row)
# #                 )

# #             rates[demographics[row[cols.index("seqn")]]["age"]][overall] += 1

# # for age in sorted(rates):
# #     positives = rates[age]["1"]
# #     negatives = rates[age]["2"]
# #     total = positives + negatives
# #     print("%s\t%s\t%s\t%s" % (age, positives, negatives, total))
