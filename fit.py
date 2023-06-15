#!/usr/bin/env python3
from pathlib import Path
from textwrap import dedent

import numpy as np
import pandas as pd

import stats
from mgs import MGSData, TaxID, rna_bioprojects
from pathogens import iter_pathogens


def geom_mean(x: pd.DataFrame) -> float:
    return np.exp(np.mean(np.log(x)))


def print_summary(
    pathogen: str,
    taxids: frozenset[TaxID],
    predictor: str,
    study: str,
    ra_per_predictor: pd.DataFrame,
) -> None:
    title = f"{pathogen} {list(taxids)} relative abundance per {predictor} in {study.title()}"
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
    for pathogen_name, predictor_type, taxids, predictors in iter_pathogens():
        for study, bioproject in rna_bioprojects.items():
            model = stats.build_model(
                mgs_data,
                bioproject,
                predictors,
                taxids,
                random_seed=1,
            )
            model.fit_model()
            df = model.dataframe
            assert df is not None
            ra_per_predictor = pd.pivot_table(
                df, index="draws", values=["ra_per_predictor"]
            )
            print_summary(
                pathogen_name, taxids, predictor_type, study, ra_per_predictor
            )
            continue
            fig_hist = model.plot_posterior_histograms()
            fig_hist.savefig(figdir / f"{study}-{pathogen_name}-posthist.pdf")
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
                g.savefig(figdir / f"{study}-{pathogen_name}-{y}-vs-{x}.pdf")

            df.to_csv(
                outdir / f"{study}-{pathogen_name}.tsv.gz",
                sep="\t",
                index=False,
                compression="gzip",
            )


if __name__ == "__main__":
    start()
