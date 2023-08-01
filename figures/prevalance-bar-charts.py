#!/usr/bin/env python3

import sys

sys.path.append("..")

from collections import defaultdict

import matplotlib.pyplot as plt  # type: ignore
import matplotlib.ticker as mtick  # type: ignore
import numpy as np

import pathogens


def start():
    labels = []
    us = []
    dk = []
    scores = []

    for (
        pathogen_name,
        tidy_name,
        predictor_type,
        taxids,
        predictors,
    ) in pathogens.predictors_by_taxid():
        (taxid,) = taxids

        if predictor_type != "prevalence":
            continue

        us_predictor = None
        dk_predictor = None
        for predictor in predictors:
            if predictor.state:
                continue

            if predictor.country == "Denmark":
                if 2015 <= predictor.parsed_start.year <= 2018:
                    if dk_predictor is None:
                        dk_predictor = predictor
                    elif (
                        dk_predictor.infections_per_100k
                        != predictor.infections_per_100k
                    ):
                        raise Exception("%s vs %s" % (dk_predictor, predictor))
            elif predictor.country == "United States":
                if 2020 <= predictor.parsed_start.year <= 2021:
                    if us_predictor is None:
                        us_predictor = predictor
                    elif (
                        us_predictor.infections_per_100k
                        != predictor.infections_per_100k
                    ):
                        raise Exception("%s vs %s" % (us_predictor, predictor))

        assert us_predictor is not None
        assert dk_predictor is not None

        labels.append(tidy_name)
        # Convert per 100k to percent by dividing by 1000.
        us.append(us_predictor.infections_per_100k / 1000)
        dk.append(dk_predictor.infections_per_100k / 1000)

        scores.append(
            us_predictor.infections_per_100k + dk_predictor.infections_per_100k
        )

    # Sort rows by average prevalence.
    labels = [x for _, x in sorted(zip(scores, labels))]
    us = [x for _, x in sorted(zip(scores, us))]
    dk = [x for _, x in sorted(zip(scores, dk))]

    fig, ax = plt.subplots(constrained_layout=True)
    plt.title("Estimated prevalence of persistent viral infections")

    ax.set_xlim(xmax=100)

    width = 0.4
    gap = 0.04

    y_pos_us = np.arange(len(labels))
    y_pos_dk = y_pos_us + width

    ax.barh(
        y_pos_us,
        us,
        width - gap,
        label="United States, 2020-2021\n(Crits-Christoph, Rothman, Spurbeck)",
    )
    ax.barh(y_pos_dk, dk, width - gap, label="Denmark, 2015-2018 (Brinch)")
    ax.set_yticks((y_pos_us + y_pos_dk) / 2, labels=labels, fontsize=9)

    def pretty_percentage(v):
        if v < 1:
            return "%.2f%%" % v
        if v < 10:
            return "%.1f%%" % v
        return "%.0f%%" % v

    for y_pos, x_pos in zip(
        np.concatenate((y_pos_us, y_pos_dk)) - 0.15, us + dk
    ):
        ax.text(x_pos + 0.3, y_pos, pretty_percentage(x_pos))

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    ax.legend(loc="lower right")
    ax.xaxis.set_major_formatter(mtick.PercentFormatter())

    fig.savefig("prevalance-bar-charts.png", dpi=300)
    plt.clf()


if __name__ == "__main__":
    start()
