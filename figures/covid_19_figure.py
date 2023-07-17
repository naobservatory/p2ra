#!/usr/bin/env python3

import sys

sys.path.append("..")

from collections import defaultdict
import matplotlib as mpl
import matplotlib.pyplot as plt  # type: ignore
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
import matplotlib.transforms as transforms


import datetime as dt
import numpy as np

import mgs
import pathogens
import stats

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42


def start():
    mgs_data = mgs.MGSData.from_repo()
    for (
        pathogen_name,
        predictor_type,
        _,
        taxids,
        predictors,
    ) in sorted(pathogens.predictors_by_taxid()):
        if pathogen_name != "sars_cov_2":
            continue
        for study, bioproject in sorted(mgs.target_bioprojects.items()):
            if study == "brinch":
                continue

            matching_reads = mgs_data.viral_reads(*bioproject, taxids)

            viral_samples = mgs_data.sample_attributes(
                *bioproject,
                enrichment=mgs.Enrichment.VIRAL,
            )
            sample_locations = set()
            sample_dates = set()

            for sample in viral_samples.values():
                for predictor in stats.lookup_variables(sample, predictors):
                    sample_dates.add(predictor.get_date())
                    sample_locations.add(
                        (sample.country, sample.state, sample.county)
                    )

            extra_dates = dt.timedelta(days=14)

            start_date = min(sample_dates) - extra_dates
            end_date = max(sample_dates) + extra_dates

            chosen_predictors = []
            for predictor in pathogens.pathogens[
                pathogen_name
            ].estimate_incidences():
                if (
                    predictor.country,
                    predictor.state,
                    predictor.county,
                ) in sample_locations:
                    if start_date <= predictor.get_date() <= end_date:
                        chosen_predictors.append(predictor)

            fig, axs = plt.subplots(
                2,
                1,
                constrained_layout=True,
                sharex=True,
                gridspec_kw={"height_ratios": [7, 1]},
            )

            for county in sorted(
                set(predictor.county for predictor in chosen_predictors)
            ):
                county_predictors = []
                for predictor in chosen_predictors:
                    if predictor.county == county:
                        county_predictors.append(
                            (predictor.get_date(), predictor.get_data() / 52)
                        )

                county_predictors = sorted(
                    county_predictors, key=lambda x: x[0]
                )
                if study == "spurbeck":
                    axs[0].plot(
                        *zip(*county_predictors),
                        # no labels
                        label="_nolegend_",
                        color="#84b3da",
                        linewidth=0.5,
                    )

                else:
                    axs[0].plot(
                        *zip(*county_predictors),
                        label=county,
                        linewidth=2,
                    )
            author_name = study.replace("_", "-").title()
            axs[0].set_title(
                f"SARS-CoV-2 Cases per Week in {author_name} et al. counties"
            )

            if study == "spurbeck":
                date_incidence_dict = defaultdict(list)
                for predictor in chosen_predictors:
                    date_incidence_dict[predictor.get_date()].append(
                        predictor.get_data() / 52
                    )

                average_incidences = [
                    (date, sum(incidences) / len(incidences))
                    for date, incidences in date_incidence_dict.items()
                ]

                average_incidences = sorted(
                    average_incidences, key=lambda x: x[0]
                )
                # drop all current labels in the legend

                axs[0].plot(
                    *zip(*average_incidences),
                    label="Average",
                    color="#2074b4",
                    linewidth=2,
                )
            axs[0].set_ylabel("Incidence per 100k per week")

            axs[0].legend(frameon=False)

            axs[1].scatter(
                sorted(list(sample_dates)),
                [0.5] * len(sample_dates),
                color="darkblue",
                alpha=0.2,
                edgecolor="none",
            )
            dates = [predictor.get_date() for predictor in chosen_predictors]
            earliest_date = min(dates)
            axs[0].set_xlim(left=earliest_date)

            axs[0].spines["right"].set_visible(False)
            axs[0].spines["top"].set_visible(False)
            axs[0].spines["bottom"].set_visible(False)
            axs[0].spines["left"].set_visible(False)

            axs[0].xaxis.set_visible(True)
            axs[0].xaxis.set_ticks_position("none")
            axs[0].axhline(0, color="black", linewidth=0.8)

            interval = (
                40
                if study == "spurbeck"
                else 10
                if study == "crits_christoph"
                else 20
            )

            ymin, ymax = axs[0].get_ylim()

            axs[0].set_yticks(np.arange(0, ymax, interval))

            for i in range(int(interval), int(ymax), int(interval)):
                axs[0].axhline(
                    i,
                    color="lightgrey",
                    linewidth=0.7,
                    linestyle="--",
                    zorder=1,
                )

            axs[0].yaxis.set_ticks_position("none")

            axs[0].set_ylabel("")

            axs[1].spines["right"].set_visible(False)
            axs[1].spines["top"].set_visible(False)
            axs[1].set_yticks([])
            axs[1].set_yticklabels([])

            axs[1].set_ylabel("")
            if study == "crits_christoph":
                axs[0].xaxis.set_major_locator(mdates.MonthLocator())
                axs[0].xaxis.set_major_formatter(
                    mdates.DateFormatter("%b\n%Y")
                )
                axs[0].xaxis.set_minor_locator(
                    mdates.WeekdayLocator(byweekday=mdates.MO)
                )
                axs[0].xaxis.set_minor_formatter(mdates.DateFormatter("%U"))
            else:
                axs[0].xaxis.set_major_locator(mdates.YearLocator())

                axs[0].xaxis.set_major_formatter(
                    mdates.DateFormatter("%b\n%Y")
                )

                axs[0].xaxis.set_minor_locator(mdates.MonthLocator())

                axs[0].xaxis.set_minor_formatter(mdates.DateFormatter("%b"))

            plt.savefig(f"{study}_{pathogen_name}_graph.pdf", dpi=300)

            plt.close(fig)


if __name__ == "__main__":
    start()
