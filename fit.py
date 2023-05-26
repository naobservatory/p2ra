#!/usr/bin/env python3
from pathlib import Path
from textwrap import dedent

import numpy as np
import pandas as pd

import stats
from mgs import BioProject, MGSData


def geom_mean(x: np.ndarray) -> float:
    return np.exp(np.mean(np.log(x)))


def print_summary(
    study: str, pathogen: str, predictor: str, ra_per_predictor: np.ndarray
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


bioprojects = {
    "crits_christoph": BioProject("PRJNA661613"),
    "rothman": BioProject("PRJNA729801"),
}


def start():
    outdir = Path("fits")
    outdir.mkdir(exist_ok=True)
    mgs_data = MGSData.from_repo()
    predictor = "incidence"
    for pathogen_name in ["sars_cov_2", "norovirus"]:
        for study, bioproject in bioprojects.items():
            model = stats.build_model(
                mgs_data, bioproject, pathogen_name, predictor
            )
            model.fit_model(random_seed=1)
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
            print_summary(study, pathogen_name, predictor, ra_per_predictor)


if __name__ == "__main__":
    start()