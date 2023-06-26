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
        predictor_type,
        taxids,
        predictors,
    ) in pathogens.predictors_by_taxid():
        pathogen_taxids_by_name[pathogen_name].append(taxids)

    row_names = []
    row_scores = []
    col_names = [
        "-".join(word.capitalize() for word in author.split("_"))
        for author in sorted(mgs.rna_bioprojects)
    ]
    # columns: studies
    # rows: pathogens
    table_text = []
    table_colors = []

    n_colors = 6
    colors = plt.cm.BuPu(np.linspace(0, 0.4, n_colors))

    for (
        pathogen_name,
        predictor_type,
        taxids,
        predictors,
    ) in sorted(pathogens.predictors_by_taxid()):
        name = tidy_name(pathogen_name, taxids)
        row_names.append(name)

        table_row_text = []
        table_row_colors = []

        row_score = 0

        for study, bioproject in sorted(mgs.rna_bioprojects.items()):
            matching_reads = mgs_data.viral_reads(bioproject, taxids)
            viral_samples = mgs_data.sample_attributes(
                bioproject,
                enrichment=mgs.Enrichment.VIRAL,
            )

            n_samples = len(viral_samples)
            counts = [
                count
                for sample, count in matching_reads.items()
                if sample in viral_samples and count
            ]
            n_samples_with_match = len(counts)
            n_matches = sum(counts)

            ratio = n_samples_with_match / n_samples
            if n_samples_with_match == 0:
                color = (1, 1, 1)
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
        row_scores.append((-row_score, name))

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

    fig, ax = plt.subplots(constrained_layout=True, figsize=(6, 3))
    ax.set_axis_off()

    table = ax.table(
        cellText=table_text,
        cellColours=table_colors,
        rowLabels=row_names,
        colLabels=col_names,
        cellLoc="center",
        loc="upper left",
    )
    fig.savefig("matching-reads-table.png", dpi=180)
    plt.clf()


if __name__ == "__main__":
    start()
