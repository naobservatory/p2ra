#!/usr/bin/env python3
from pathlib import Path

import pandas as pd
import numpy as np
import numpy.typing as npt

import stats
from mgs import MGSData, rna_bioprojects
from pathogens import predictors_by_taxid

PERCENTILES = [5, 25, 50, 75, 95]


def geom_mean(x: npt.ArrayLike) -> float:
    return np.exp(np.mean(np.log(x)))


# TODO: Deprecate
def summarize_output(coeffs: pd.DataFrame) -> pd.DataFrame:
    # for location in self.fine_locations + ["Overall"]:
    #     c = coeffs[location]
    #     output.append(
    #         [location, np.mean(c), geom_mean(c)]
    #         + list(np.percentile(c, PERCENTILES))
    #     )
    summary = pd.DataFrame()
    return summary


def start() -> None:
    figdir = Path("fig")
    figdir.mkdir(exist_ok=True)
    mgs_data = MGSData.from_repo()
    output_data = []
    for (
        pathogen_name,
        predictor_type,
        taxids,
        predictors,
    ) in predictors_by_taxid():
        for study, bioproject in rna_bioprojects.items():
            model = stats.build_model(
                mgs_data,
                bioproject,
                predictors,
                taxids,
                random_seed=1,
            )
            model.fit_model()
            model.plot_figures(
                path=figdir,
                prefix=f"{pathogen_name}-{list(taxids)}-{predictor_type}-{study}",
            )
            out = model.get_coefficients()
            out["pathogen"] = pathogen_name
            out["taxids"] = ",".join(str(t) for t in taxids)
            out["predictor_type"] = predictor_type
            out["study"] = study
            output_data.append(out)
    coeffs = pd.concat(output_data)
    coeffs.to_csv("fits.tsv", sep="\t", index=False)
    # TODO: summarize into a separate csv?
    summary = summarize_output(coeffs)
    summary.to_csv("fits_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    start()
