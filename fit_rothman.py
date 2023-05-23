#!/usr/bin/env python3
from textwrap import dedent

import numpy as np
import pandas as pd

import stats
from mgs import BioProject, MGSData
from pathogens import pathogens


def geom_mean(x: np.ndarray) -> float:
    return np.exp(np.mean(np.log(x)))


def print_summary(pathogen: str, model_ra_per100: np.ndarray) -> None:
    title = f"{pathogen} relative abundance at 1% prevalence"
    percentiles = [5, 25, 50, 75, 95]
    percentile_values = np.percentile(model_ra_per100, percentiles)
    d = 1
    sep = " " * 4
    output = f"""
    {"-" * len(title)}
    {title}
    {"-" * len(title)}
    Posterior arithmetic mean:
    {np.mean(model_ra_per100):.{d}e}
    Posterior geometric mean:
    {geom_mean(model_ra_per100):.{d}e}
    Posterior quantiles:
    {sep.join(f"{p:>{d+5}}%" for p in percentiles)}
    {sep.join(f"{x:.{d}e}" for x in percentile_values)}
    """
    print(dedent(output))


def start():
    bioproject = BioProject("PRJNA729801")  # Rothman
    mgs_data = MGSData.from_repo()
    for pathogen_name in ["sars_cov_2", "norovirus"]:
        pathogen = pathogens[pathogen_name]
        model_df = stats.fit_model(
            mgs_data, bioproject, pathogen, random_seed=1
        )
        model_df.to_csv(
            f"fits/rothman-{pathogen_name}.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )
        model_ra_per100 = pd.pivot_table(
            model_df, index="draws", values=["ra_per_one_percent"]
        )
        print_summary(pathogen_name, model_ra_per100)


if __name__ == "__main__":
    start()
