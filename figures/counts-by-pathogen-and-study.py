#!/usr/bin/env python3

import sys

sys.path.append("..")

from collections import defaultdict

import matplotlib.pyplot as plt  # type: ignore
import numpy as np

import mgs
import pathogens


def start():
    mgs_data = mgs.MGSData.from_repo()

    pathogen_taxids_by_name = defaultdict(list)
    for (
        pathogen_name,
        _,
        _,
        taxids,
        _,
    ) in pathogens.predictors_by_taxid():
        pathogen_taxids_by_name[pathogen_name].append(taxids)

    row_names = []
    row_scores = []
    col_names = ["Nucleic Acid", "Selection"] + [
        "-".join(word.capitalize() for word in author.split("_"))
        for author in sorted(mgs.target_bioprojects)
    ]
    white_rgb = 1, 1, 1
    yellow_rgb = 1, 1, 0.94
    purple_rgb = 1, 0.94, 1
    green_rgb = 0.93, 1, 0.93
    red_rgb = 1, 0.96, 0.96

    col_colors = [white_rgb, white_rgb] + [
        yellow_rgb if author == "brinch" else purple_rgb
        for author in sorted(mgs.target_bioprojects)
    ]

    # columns: studies
    # rows: pathogens
    table_text = []
    table_colors = []

    n_colors = 6
    colors = plt.cm.BuPu(np.linspace(0, 0.4, n_colors))

    for (
        pathogen_name,
        tidy_name,
        predictor_type,
        taxids,
        _,
    ) in sorted(pathogens.predictors_by_taxid()):
        name = tidy_name
        row_names.append(name)

        white_rgb = 1, 1, 1

        na_type = pathogens.pathogens[
            pathogen_name
        ].pathogen_chars.na_type.value
        selection = pathogens.pathogens[
            pathogen_name
        ].pathogen_chars.selection.value
        table_row_text = [na_type, selection]

        table_row_colors = [
            yellow_rgb if na_type == "DNA" else purple_rgb,
            green_rgb if selection == "Round 1" else red_rgb,
        ]

        row_score = 0

        for study, bioprojects in sorted(mgs.target_bioprojects.items()):
            n_samples = 0
            n_samples_with_match = 0
            n_matches = 0
            for bioproject in bioprojects:
                matching_reads = mgs_data.viral_reads(bioproject, taxids)
                viral_samples = mgs_data.sample_attributes(
                    bioproject,
                    enrichment=None
                    if study == "brinch"
                    else mgs.Enrichment.VIRAL,
                )
                n_samples += len(viral_samples)
                counts = [
                    count
                    for sample, count in matching_reads.items()
                    if sample in viral_samples and count
                ]
                n_samples_with_match += len(counts)
                n_matches += sum(counts)

            ratio = n_samples_with_match / n_samples
            if n_samples_with_match == 0:
                color = white_rgb
            elif n_samples_with_match == 1:
                color = colors[1]
            elif ratio < 1 / 55:
                color = colors[2]
            elif ratio < 1 / 7:
                color = colors[3]
            elif ratio < 1 / 2:
                color = colors[4]
            else:
                color = colors[5]

            table_row_colors.append(tuple(color))
            table_row_text.append(
                "%s/%s (%s)" % (n_samples_with_match, n_samples, n_matches)
            )

            row_score += n_samples_with_match

        table_text.append(table_row_text)
        table_colors.append(table_row_colors)
        # include name as a tiebreaker, for consistent sorting below
        row_scores.append(
            (1 if na_type == "DNA" else 0, selection, -row_score, name)
        )

    # now sort everything by row_scores
    assert (
        len(row_scores)
        == len(row_names)
        == len(table_text)
        == len(table_colors)
    )
    row_names = [x for _, x in sorted(zip(row_scores, row_names))]
    table_text = [x for _, x in sorted(zip(row_scores, table_text))]
    table_colors = [x for _, x in sorted(zip(row_scores, table_colors))]

    fig, ax = plt.subplots(constrained_layout=True, figsize=(6, 3.7))
    ax.set_axis_off()

    table = ax.table(
        cellText=table_text,
        cellColours=table_colors,
        colColours=col_colors,
        rowLabels=row_names,
        colLabels=col_names,
        cellLoc="center",
        loc="upper left",
    )
    fig.savefig("matching-reads-table.png", dpi=180)
    plt.clf()


if __name__ == "__main__":
    start()
