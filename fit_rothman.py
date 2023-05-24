#!/usr/bin/env python3
from textwrap import dedent

import numpy as np
import pandas as pd

import stats
from mgs import BioProject, Enrichment, MGSData
from pathogen_properties import Predictor
from pathogens import pathogens


def geom_mean(x: np.ndarray) -> float:
    return np.exp(np.mean(np.log(x)))


def print_summary(pathogen: str, ra_per_predictor: np.ndarray) -> None:
    title = f"{pathogen} relative abundance per predictor"
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


def start():
    bioproject = BioProject("PRJNA729801")  # Rothman
    mgs_data = MGSData.from_repo()
    for pathogen_name in ["sars_cov_2", "norovirus"]:
        pathogen = pathogens[pathogen_name]
        taxids = pathogen.pathogen_chars.taxids
        samples = mgs_data.sample_attributes(
            bioproject, enrichment=Enrichment.VIRAL
        )
        incidences: list[Predictor] = pathogen.estimate_incidences()
        data = [
            stats.DataPoint(
                sample=s,
                attrs=attrs,
                viral_reads=mgs_data.viral_reads(bioproject, taxids)[s],
                predictor=stats.lookup_variable(attrs, incidences),
            )
            for s, attrs in samples.items()
        ]
        model = stats.Model(data=data)
        model.fit_model(random_seed=1)
        df = model.get_fit_dataframe()
        df.to_csv(
            f"fits/rothman-{pathogen_name}.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )
        ra_per_predictor = pd.pivot_table(
            df, index="draws", values=["ra_per_predictor"]
        )
        print_summary(pathogen_name, ra_per_predictor)


if __name__ == "__main__":
    start()
