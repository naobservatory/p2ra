# Script to create a tsv file based on the pathogen *.py files

# TODOs:
# 1. Review the column swapping logic and confirm if it meets your requirements.
# 2. Make code generalizable to other pathogen files.
# 3. Add error handling and/or logging to make debugging easier if necessary.
# 4. Add comments to explain the code.
# 5. Have it reviewed by Jeff


import csv

from hsv_1_prevalence import *

# Creating header from keys

prevalence_vars_keys = []
for outer_key in prevalence_vars.keys():
    for inner_key in prevalence_vars[outer_key].keys():
        if inner_key not in prevalence_vars_keys:
            prevalence_vars_keys.append(inner_key)

pathogen_chars_keys = list(pathogen_chars.keys())

header = ["variable_id"] + pathogen_chars_keys + prevalence_vars_keys

# Create a list to store rows
rows = []

rows.append(header)

# Writing rows
for var_id, var_data in prevalence_vars.items():
    row = [var_id]
    for key in header[1:]:
        if key in pathogen_chars:
            row.append(pathogen_chars[key])
        else:
            row.append(var_data.get(key, ""))
    rows.append(row)

# Swapping columns
variable_id_index = 0
pathogen_name_index = 1
for row in rows:
    row[variable_id_index], row[pathogen_name_index] = (
        row[pathogen_name_index],
        row[variable_id_index],
    )

# Write rows to the file
with open("prevalence_vars.tsv", "w", newline="", encoding="utf-8") as tsvfile:
    writer = csv.writer(tsvfile, delimiter="\t")
    writer.writerows(rows)
