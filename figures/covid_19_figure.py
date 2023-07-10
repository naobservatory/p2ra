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

        for study, bioproject in sorted(mgs.load_bioprojects.items()):
            matching_reads = mgs_data.viral_reads(bioproject, taxids)
            viral_samples = mgs_data.sample_attributes(
                bioproject,
                enrichment=mgs.Enrichment.VIRAL,
            )
            chosen_predictors = stats.lookup_variables(
                viral_samples, predictors
            )
            fig, ax = plt.subplots(constrained_layout=True)

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

                axes = ax.plot(*zip(*county_predictors), label=county)
                color = axes[0].get_color()

            ax.set_title(f"{pathogen_name} Predictors by County for {study}")
            ax.set_xlabel("Date")
            ax.set_ylabel("Predictor data")
            ax.legend()

            plt.savefig(f"{bioproject}_{pathogen_name}_graph.png")
            print("this worked")
            plt.close(fig)

        #     ax

        # incidences_by_day = []
        # for incidence in pathogens["sars_cov_2"].estimate_incidences():
        #     if not incidence.county:
        #         incidences_by_day.append(
        #             (
        #                 incidence.get_date(),
        #                 incidence.annual_infections_per_100k / 52,
        #             )
        #         )


if __name__ == "__main__":
    start()
