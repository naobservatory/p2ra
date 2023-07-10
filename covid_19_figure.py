#!/usr/bin/env python3

import sys

sys.path.append("..")

from collections import defaultdict

import matplotlib.pyplot as plt  # type: ignore
import numpy as np

import mgs
import pathogens
import stats


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
            print(bioproject)
            print(type(bioproject))
            matching_reads = mgs_data.viral_reads(
                bioproject, taxids
            )  # Have to account for bioproject being a list now. Discuss with Jeff tomorrow.

            viral_samples = mgs_data.sample_attributes(
                bioproject,
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
                gridspec_kw={"height_ratios": [4, 1]},
            )

            for county in sorted(
                set(predictor.county for predictor in chosen_predictors)
            ):
                county_predictors = []
                for predictor in chosen_predictors:
                    if predictor.county == county:
                        county_predictors.append(
                            (predictor.get_date(), predictor.get_data())
                        )

                county_predictors = sorted(
                    county_predictors, key=lambda x: x[0]
                )

                axs[0].plot(*zip(*county_predictors), label=county)

            axs[0].set_title(
                f"SARS-CoV-2 Incidence in {study.capitalize()} et al."
            )

            axs[0].set_ylabel("Incidence per 100k per week")

            axs[0].spines["right"].set_visible(False)
            axs[0].spines["top"].set_visible(False)

            if study == "crits_christoph":
                axs[0].xaxis.set_major_locator(
                    mdates.DayLocator(bymonthday=(1, 7, 14, 21, 28))
                )
                axs[0].xaxis.set_major_formatter(
                    mdates.DateFormatter("%Y-%m-%d")
                )
            else:
                axs[0].xaxis.set_major_locator(mdates.MonthLocator())
                axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))

            axs[0].legend(frameon=False)

            axs[1].scatter(
                sorted(list(sample_dates)),
                # set y coordiante to 1.5 for all points
                [1.5] * len(sample_dates),
                color="darkred",
                marker="|",
            )  # Here you need to customize this plot according to your requirements

            axs[1].spines["right"].set_visible(False)
            axs[1].spines["top"].set_visible(False)

            axs[0].spines["right"].set_visible(False)
            axs[0].spines["top"].set_visible(False)
            axs[0].spines["bottom"].set_visible(False)
            # remove ticks on the y axis
            axs[1].set_yticks([])
            axs[1].set_yticklabels([])
            axs[1].set_ylabel("Sample\ndates")
            # hide x axis of subplot 0
            axs[0].xaxis.set_visible(False)
            # line at y = 0
            axs[0].axhline(0, color="black", linewidth=0.8)

            plt.xticks(rotation=45)
            plt.savefig(f"{study}_{pathogen_name}_graph.png", dpi=300)

            plt.close(fig)


if __name__ == "__main__":
    start()
