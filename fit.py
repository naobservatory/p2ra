#!/usr/bin/env python3
from pathlib import Path

import pandas as pd

import stats
from mgs import Enrichment, MGSData, target_bioprojects
from pathogens import predictors_by_taxid


def summarize_output(coeffs: pd.DataFrame) -> pd.DataFrame:
    return coeffs.groupby(
        ["pathogen", "taxids", "predictor_type", "study", "location"]
    ).ra_at_1in1000.describe(percentiles=[0.05, 0.25, 0.5, 0.75, 0.95])


def start(num_samples: int, plot: bool) -> None:
    figdir = Path("fig")
    if plot:
        figdir.mkdir(exist_ok=True)
    mgs_data = MGSData.from_repo()
    input_data = []
    output_data = []
    for (
        _,
        pathogen_name,
        predictor_type,
        taxids,
        predictors,
    ) in predictors_by_taxid():
        taxids_str = "_".join(str(t) for t in taxids)
        for study, bioprojects in target_bioprojects.items():
            enrichment = None if study == "brinch" else Enrichment.VIRAL
            model = stats.build_model(
                mgs_data,
                bioprojects,
                predictors,
                taxids,
                random_seed=sum(taxids),
                enrichment=enrichment,
            )
            model.fit_model(num_samples=num_samples)
            if plot:
                model.plot_figures(
                    path=figdir,
                    prefix=f"{pathogen_name}-{predictor_type}-{study}",
                )
            metadata = dict(
                pathogen=pathogen_name,
                taxids=taxids_str,
                predictor_type=predictor_type,
                study=study,
            )
            input_data.append(model.input_df.assign(**metadata))
            output_data.append(model.get_coefficients().assign(**metadata))
    input = pd.concat(input_data)
    input.to_csv("input.tsv", sep="\t", index=False)
    coeffs = pd.concat(output_data)
    coeffs.to_csv("fits.tsv", sep="\t", index=False)
    summary = summarize_output(coeffs)
    summary.to_csv("fits_summary.tsv", sep="\t")


if __name__ == "__main__":
    # TODO: Command line arguments
    start(num_samples=1000, plot=True)
