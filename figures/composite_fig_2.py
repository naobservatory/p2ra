import sys as sys
import csv

sys.path.append("..")

import matplotlib as mpl
import matplotlib.pyplot as plt  # type: ignore
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
from matplotlib.gridspec import GridSpec


import datetime as dt
import numpy as np
import pandas as pd
from dataclasses import dataclass


import mgs
import pathogens

from populations import us_population

sys.path.append("..")

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42


def transform_study_name(s):
    return "-".join([word.capitalize() for word in s.split("_")])


def get_averaged_incidence(df, pathogen_name):
    averaged_incidence = []

    for date, group in df[df["pathogen"] == pathogen_name].groupby("date"):
        weights = group.apply(
            lambda row: us_population(
                row["date"].year,
                row["county"],
                row["state"],
            ).people,
            axis=1,
        )

        incidence = np.average(group["incidence"], weights=weights)
        averaged_incidence.append((date, incidence))

    return pd.DataFrame(
        averaged_incidence, columns=["date", f"{pathogen_name}"]
    )


def assemble_incidence_data():
    target_pathogens = ["sars_cov_2", "influenza", "norovirus"]

    norovirus_added_days = dt.timedelta(days=15.5)
    influenza_added_days = dt.timedelta(days=3.5)

    pathogen_predictors = []
    for pathogen_name in target_pathogens:
        for predictor in pathogens.pathogens[
            pathogen_name
        ].estimate_incidences():
            if pathogen_name == "norovirus":
                date = predictor.get_dates()[0] + norovirus_added_days

            elif pathogen_name == "influenza":
                date = predictor.get_dates()[0] + influenza_added_days

            elif pathogen_name == "sars_cov_2":
                date = predictor.get_dates()[
                    0
                ]  # daily, so no need to add days to get to the middle of the time period

            location = (*predictor.get_location(),)  # country, state, county
            incidence = predictor.get_data()

            if date < dt.date(2020, 5, 15) or date > dt.date(2022, 1, 15):
                continue

            pathogen_predictors.append(
                (
                    pathogen_name,
                    date,
                    *location,
                    incidence,
                )
            )

    df = pd.DataFrame(
        pathogen_predictors,
        columns=[
            "pathogen",
            "date",
            "country",
            "state",
            "county",
            "incidence",
        ],
    )

    # Norovirus and Influenza have subtypes that we want to aggregate on, hence we are summing estimates with the same pathogen, date and location.
    df = df.groupby(
        ["pathogen", "country", "state", "county", "date"],
        as_index=False,
        dropna=False,
    ).agg({"incidence": "sum"})

    # sort by date
    df = df.sort_values(by=["date"]).fillna("")

    df_final = (
        get_averaged_incidence(df, "sars_cov_2")
        .merge(
            get_averaged_incidence(df, "influenza"),
            on="date",
            how="outer",
        )
        .merge(
            # no need to average, as data is already across the US
            df[df["pathogen"] == "norovirus"][["date", "incidence"]].rename(
                columns={"incidence": "norovirus"}
            ),
            on="date",
            how="outer",
        )
    )
    return df_final


def prevalence_bar_chart(ax1, barplot_colors):
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
    ax1.set_xlim(xmax=100)

    width = 0.4
    gap = 0.04

    y_pos_dk = np.arange(len(labels))
    y_pos_us = y_pos_dk + width

    ax1.barh(
        y_pos_us,
        us,
        width - gap,
        color=barplot_colors[0],
        label="United States, 2020-2021\n(Crits-Christoph, Rothman, Spurbeck)",
    )
    ax1.barh(
        y_pos_dk,
        dk,
        width - gap,
        label="Denmark, 2015-2018 (Brinch)",
        color=barplot_colors[1],
    )
    ax1.set_yticks((y_pos_us + y_pos_dk) / 2, labels=labels, fontsize=9)

    def pretty_percentage(v):
        if v < 1:
            return "%.2f%%" % v
        if v < 10:
            return "%.1f%%" % v
        return "%.0f%%" % v

    for y_pos, x_pos in zip(
        np.concatenate((y_pos_us, y_pos_dk)) - 0.15, us + dk
    ):
        ax1.text(x_pos + 0.3, y_pos, pretty_percentage(x_pos))

    ax1.spines["right"].set_visible(False)
    ax1.spines["top"].set_visible(False)

    ax1.legend(loc="lower right")
    ax1.xaxis.set_major_formatter(mticker.PercentFormatter())
    ax1_title = ax1.set_title(
        "a",
        fontweight="bold",
        loc="left",
    )

    ax1_title.set_position((-0.08, 0))

    return ax1


def get_sample_dates():
    mgs_data = mgs.MGSData.from_repo()
    sample_dates = {
        "rothman": [],
        "crits_christoph": [],
        "spurbeck": [],
    }
    for study, bioproject in sorted(mgs.target_bioprojects.items()):
        if study in ["brinch"]:
            continue

        viral_samples = mgs_data.sample_attributes(
            *bioproject,
            enrichment=mgs.Enrichment.VIRAL,
        )
        sum = 0
        for sample, sample_attributes in viral_samples.items():
            sample_dates[study].append(sample_attributes.date)
    return sample_dates


def all_incidence_figure(
    df,
    sample_dates,
    ax2,
    ax3,
    transform_study_name,
    virus_plot_colors,
    study_coverage_colors,
):
    figure_labels = {
        "sars_cov_2": "SARS-CoV-2",
        "influenza": "Influenza",
        "norovirus": "Norovirus",
    }

    ax2.plot(
        df["date"],
        df["sars_cov_2"],
        label=figure_labels["sars_cov_2"],
        color=virus_plot_colors["sars_cov_2"],
    )

    # ax2.plot(
    #     df.dropna(subset=["influenza"]).sort_values(by="date")["date"],
    #     df.dropna(subset=["influenza"]).sort_values(by="date")["influenza"],
    #     label=figure_labels["influenza"],
    #     color=virus_plot_colors["influenza"],
    # )

    # ax2.plot(
    #     df.dropna(subset=["norovirus"]).sort_values(by="date")["date"],
    #     df.dropna(subset=["norovirus"]).sort_values(by="date")["norovirus"],
    #     label=figure_labels["norovirus"],
    #     color=virus_plot_colors["norovirus"],
    # )

    study_coverage_y_values = range(1, len(sample_dates) + 1)

    for study, dates in sample_dates.items():
        ax2.axvspan(
            min(dates),
            max(dates) + dt.timedelta(days=2),
            ymin=0.035,
            facecolor=study_coverage_colors[study],
            alpha=0.1,
        )

    for index, (study, dates) in zip(
        study_coverage_y_values, sample_dates.items()
    ):
        ax3.plot(
            [
                min(dates),
                max(dates) + dt.timedelta(days=2),
            ],  # hacky, have to add two days so Crits-Cristoph shows up in the study coverage plot
            [index, index],
            color=study_coverage_colors[study],
            # alpha = 0.4,
            linewidth=5,
            solid_capstyle="round",
        )

        x_offset = dt.timedelta(days=8)
        ax3.text(
            max(dates) + x_offset,
            index,
            transform_study_name(study),
            fontsize=10.5,
            va="center",
            # color=study_coverage_colors[study],
        )
    ax2.spines["right"].set_visible(False)
    ax2.spines["top"].set_visible(False)
    ax2.spines["bottom"].set_visible(False)
    ax2.spines["left"].set_visible(False)

    ax2.xaxis.set_visible(True)
    ax2.xaxis.set_ticks_position("none")
    plt.setp(ax2.get_xticklabels(), visible=False)
    ax2.axhline(0, color="black", linewidth=0.8)

    ax2.yaxis.set_ticks_position("none")

    ax2_title = ax2.set_title(
        "b",
        fontweight="bold",
        loc="left",
    )

    ax2_title.set_position((-0.08, 0))
    x_min, _ = ax2.get_xlim()
    last_date = df.dropna(
        subset=["sars_cov_2", "influenza", "norovirus"], how="all"
    )["date"].max()

    x_axis_max_offset = dt.timedelta(days=7)

    ax2.set_xlim(x_min, last_date + x_axis_max_offset)

    for pathogen in ["sars_cov_2", "influenza", "norovirus"]:
        last_incidence_value = (
            df.sort_values(by="date")
            .dropna(subset=[pathogen])[pathogen]
            .iloc[-1]
        )

        last_date = ax2.get_xlim()[1]

        ax2.text(
            last_date,
            last_incidence_value,
            figure_labels[pathogen],
            color=virus_plot_colors[pathogen],
            verticalalignment="center",
            fontsize=10,
        )

    ax3.text(  # title
        dt.datetime(2020, 4, 16),
        5,
        "Study coverage",
        fontsize=12,
        va="center",
    )

    ax2.yaxis.set_major_formatter(
        mticker.FuncFormatter(lambda x, p: format(int(x), ","))
    )

    ax3.spines["right"].set_visible(False)
    ax3.spines["top"].set_visible(False)
    ax3.set_ylim(0, max(study_coverage_y_values) + 1)
    ax3.set_yticks([])
    ax3.set_yticklabels([])

    for y_loc in range(2000, 10001, 2000):
        ax2.axhline(y_loc, color="gray", linestyle="-", lw=0.5, alpha=0.5)

    ax3.xaxis.set_major_locator(mdates.YearLocator())
    ax3.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))

    ax3.xaxis.set_minor_locator(mdates.MonthLocator())
    ax3.xaxis.set_minor_formatter("")

    # ax2.legend(loc="upper left", borderaxespad=0)

    return ax2, ax3


def norovirus_influenza_figure(
    df,
    sample_dates,
    ax4,
    ax5,
    transform_study_name,
    virus_plot_colors,
    study_coverage_colors,
):
    figure_labels = {
        "sars_cov_2": "SARS-CoV-2",
        "influenza": "Influenza",
        "norovirus": "Norovirus",
    }

    ax4.plot(
        df.dropna(subset=["influenza"]).sort_values(by="date")["date"],
        df.dropna(subset=["influenza"]).sort_values(by="date")["influenza"],
        marker="o",
        markersize=2.5,
        label=f"{figure_labels['influenza']}\n(weekly data)",
        color=virus_plot_colors["influenza"],
    )

    ax4.plot(
        df.dropna(subset=["norovirus"]).sort_values(by="date")["date"],
        df.dropna(subset=["norovirus"]).sort_values(by="date")["norovirus"],
        marker="o",
        markersize=3.5,
        label=f"{figure_labels['norovirus']}\n(weekly data)",
        color=virus_plot_colors["norovirus"],
    )

    study_coverage_y_values = range(1, len(sample_dates) + 1)

    for study, dates in sample_dates.items():
        ax4.axvspan(
            min(dates),
            max(dates) + dt.timedelta(days=2),
            ymin=0.04,
            facecolor=study_coverage_colors[study],
            alpha=0.1,
        )

    for index, (study, dates) in zip(
        study_coverage_y_values, sample_dates.items()
    ):
        ax5.plot(
            [
                min(dates),
                max(dates) + dt.timedelta(days=2),
            ],  # hacky, have to add two days so Crits-Cristoph shows up in the study coverage plot
            [index, index],
            color=study_coverage_colors[study],
            # alpha=0.2,
            linewidth=5,
            solid_capstyle="round",
        )

        x_offset = dt.timedelta(days=8)
        ax5.text(
            max(dates) + x_offset,
            index,
            transform_study_name(study),
            fontsize=10.5,
            va="center",
            # color=study_coverage_colors[study],
        )

    ax4.spines["right"].set_visible(False)
    ax4.spines["top"].set_visible(False)
    ax4.spines["bottom"].set_visible(False)
    ax4.spines["left"].set_visible(False)

    ax4.xaxis.set_visible(True)
    ax4.xaxis.set_ticks_position("none")
    plt.setp(ax4.get_xticklabels(), visible=False)
    ax4.axhline(0, color="black", linewidth=0.8)
    # remove ticks on ax4 y axis, but still show labels on y axis
    ax4.yaxis.set_ticks_position("none")

    ax4_title = ax4.set_title(
        "c",
        fontweight="bold",
        loc="left",
    )

    ax4_title.set_position((-0.08, 0))

    x_min, _ = ax4.get_xlim()
    last_date = df.dropna(subset=["influenza", "norovirus"], how="all")[
        "date"
    ].max()

    x_axis_max_offset = dt.timedelta(days=7)

    ax4.set_xlim(x_min, last_date + x_axis_max_offset)

    for pathogen in ["influenza", "norovirus"]:
        last_incidence_value = (
            df.sort_values(by="date")
            .dropna(subset=[pathogen])[pathogen]
            .iloc[-1]
        )

        last_date = ax4.get_xlim()[1]

        ax4.text(
            last_date,
            last_incidence_value,
            figure_labels[pathogen],
            color=virus_plot_colors[pathogen],
            verticalalignment="center",
            fontsize=10,
        )

    ax5.text(
        dt.datetime(2020, 4, 16),
        5,
        "Study coverage",
        fontsize=12,
        va="center",
    )

    ax5.spines["right"].set_visible(False)
    ax5.spines["top"].set_visible(False)
    ax5.set_ylim(0, max(study_coverage_y_values) + 1)
    ax5.set_yticks([])
    ax5.set_yticklabels([])

    for y_loc in range(20, 161, 20):
        ax4.axhline(y_loc, color="gray", linestyle="-", lw=0.5, alpha=0.5)

    ax5.xaxis.set_major_locator(mdates.YearLocator())
    ax5.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax5.xaxis.set_minor_locator(mdates.MonthLocator())
    ax5.xaxis.set_minor_formatter("")

    # ax4.legend(loc="upper left", borderaxespad=0)

    return ax4, ax5


def start():
    df = assemble_incidence_data()
    sample_dates = get_sample_dates()

    virus_plot_colors = {
        "sars_cov_2": "#b13607",
        "influenza": "#cf0a66",
        "norovirus": "#2c8565",
    }

    study_coverage_colors = {
        "crits_christoph": "#2ca02c",
        "rothman": "#ff7f0e",
        "spurbeck": "#1f77b4",
    }

    barplot_colors = ["#8dd3c7", "#fdb462"]

    fig = plt.figure(
        figsize=(9, 12),
    )
    gs = fig.add_gridspec(5, 2, height_ratios=[7, 5, 1, 5, 1])

    ax1 = fig.add_subplot(gs[0, :])
    ax3 = fig.add_subplot(gs[2, :])
    ax2 = fig.add_subplot(gs[1, :], sharex=ax3)
    ax5 = fig.add_subplot(gs[4, :])
    ax4 = fig.add_subplot(gs[3, :], sharex=ax5)

    ax1 = prevalence_bar_chart(ax1, barplot_colors)
    ax2, ax3 = all_incidence_figure(
        df,
        sample_dates,
        ax2,
        ax3,
        transform_study_name,
        virus_plot_colors,
        study_coverage_colors,
    )
    ax4, ax5 = norovirus_influenza_figure(
        df,
        sample_dates,
        ax4,
        ax5,
        transform_study_name,
        virus_plot_colors,
        study_coverage_colors,
    )

    plt.tight_layout()
    plt.savefig(
        "composite_fig_2.pdf",
        bbox_inches="tight",
    )
    plt.savefig(
        "composite_fig_2.png",
        bbox_inches="tight",
        dpi=600,
    )


if __name__ == "__main__":
    start()
