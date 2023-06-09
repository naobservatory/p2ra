#!/usr/bin/env python3
from pathlib import Path
from textwrap import dedent

import numpy as np
import pandas as pd

import pathogens
import stats
from mgs import MGSData, rna_bioprojects
from pathogen_properties import Predictor, by_taxids


def geom_mean(x: pd.DataFrame) -> float:
    return np.exp(np.mean(np.log(x)))


def print_summary(
    study: str, pathogen: str, predictor: str, ra_per_predictor: pd.DataFrame
) -> None:
    title = f"{study.title()}: {pathogen} relative abundance per {predictor}"
    percentiles = [5, 25, 50, 75, 95]
    percentile_values = np.percentile(ra_per_predictor, percentiles)
    d = 1
    sep = " " * 4
    output = f"""
    {"-" * len(title)}
    {title}
    {"-" * len(title)}
    Posterior arithmetic mean:
    {np.mean(ra_per_predictor):.{d}e}
    Posterior geometric mean:
    {geom_mean(ra_per_predictor):.{d}e}
    Posterior quantiles:
    {sep.join(f"{p:>{d+5}}%" for p in percentiles)}
    {sep.join(f"{x:.{d}e}" for x in percentile_values)}
    """
    print(dedent(output))


def start() -> None:
    outdir = Path("fits")
    outdir.mkdir(exist_ok=True)
    figdir = Path("fig")
    figdir.mkdir(exist_ok=True)
    mgs_data = MGSData.from_repo()
    predictor = "incidence"
    for pathogen_name in ["sars_cov_2"]:
        pathogen = pathogens.pathogens[pathogen_name]
        for study, bioproject in rna_bioprojects.items():
            predictors: list[Predictor]
            if predictor == "incidence":
                predictors = pathogen.estimate_incidences()
            elif predictor == "prevalence":
                predictors = pathogen.estimate_prevalences()
            else:
                raise ValueError(
                    f"{predictor} must be one of 'incidence' or 'prevalence'"
                )

            for taxids, grouped_predictors in by_taxids(
                pathogen.pathogen_chars, predictors
            ).items():
                model = stats.build_model(
                    mgs_data,
                    bioproject,
                    grouped_predictors,
                    taxids,
                    random_seed=1,
                )
                model.fit_model()
                fig_hist = model.plot_posterior_histograms()
                fig_hist.savefig(
                    figdir / f"{study}-{pathogen_name}-posthist.pdf"
                )
                xys = [
                    ("date", "viral_reads"),
                    ("date", "predictor"),
                    ("predictor", "viral_reads"),
                ]
                for x, y in xys:
                    g = model.plot_posterior_samples(
                        x, y, style="county", hue="fine_location"
                    )
                    if y == "predictor":
                        g.set(yscale="log")
                    g.savefig(
                        figdir / f"{study}-{pathogen_name}-{y}-vs-{x}.pdf"
                    )

                df = model.dataframe
                assert df is not None
                df.to_csv(
                    outdir / f"{study}-{pathogen_name}.tsv.gz",
                    sep="\t",
                    index=False,
                    compression="gzip",
                )
                ra_per_predictor = pd.pivot_table(
                    df, index="draws", values=["ra_per_predictor"]
                )
                print_summary(
                    study, pathogen_name, predictor, ra_per_predictor
                )


if __name__ == "__main__":
    start()
