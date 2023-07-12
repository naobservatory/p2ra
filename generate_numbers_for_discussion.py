#!/usr/bin/env python3

AGNOSTIC_FOLLOWUP = 1000  # weekly reads matching pathogen
TARGETED_FOLLOWUP = 3  # weekly reads matching pathogen

SEQUENCING_COST = 8000 / 1e9  # dollars per read

cols = None
with open("fits_summary.tsv") as inf:
    for line in inf:
        row = line.strip().split("\t")
        if not cols:
            cols = row
            continue

        if row[cols.index("location")] != "Overall":
            continue

        pathogen = row[cols.index("pathogen")]
        study = row[cols.index("study")]
        median = float(row[cols.index("50%")])
        predictor_type = row[cols.index("predictor_type")]
        is_prevalence = {
            "prevalence": True,
            "incidence": False,
        }[predictor_type]

        # 1% prevalence or 0.5% weekly incidence
        adjusted_relative_abundance = median * 10 if is_prevalence else median * 5

        print(pathogen, study)
        print(
            "  Relative abundance of 1 in %.0e"
            % (1 / adjusted_relative_abundance)
        )
        print(
            "  Weekly reads to flag for manual followup: %.0e"
            % (AGNOSTIC_FOLLOWUP * 1 / adjusted_relative_abundance)
        )
        print(
            "  Weekly agnostic cost: $%.0f"
            % (
                AGNOSTIC_FOLLOWUP
                * 1
                / adjusted_relative_abundance
                * SEQUENCING_COST
            )
        )
        print(
            "  Annual agnostic cost: $%.0f"
            % (
                AGNOSTIC_FOLLOWUP
                * 1
                / adjusted_relative_abundance
                * SEQUENCING_COST
                * 52
            )
        )
        print(
            "  Annual targeted cost: $%.0f"
            % (
                TARGETED_FOLLOWUP
                * 1
                / adjusted_relative_abundance
                * SEQUENCING_COST
                * 52
            )
        )
