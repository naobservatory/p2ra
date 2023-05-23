#!/usr/bin/env python3
from textwrap import dedent

import numpy as np
import pandas as pd

import stats
from mgs import BioProject, Enrichment, MGSData
from pathogen_properties import NAType
from pathogens import pathogens


def geom_mean(x: np.ndarray) -> float:
    return np.exp(np.mean(np.log(x)))


def print_summary(
    pathogen: str, naive_ra_per100: float, model_ra_per100: np.ndarray
) -> None:
    title = f"{pathogen} relative abundance at 1% prevalence"
    percentiles = [5, 25, 50, 75, 95]
    percentile_values = np.percentile(model_ra_per100, percentiles)
    d = 1
    sep = " " * 4
    output = f"""
    {"-" * len(title)}
    {title}
    {"-" * len(title)}
    Naive estimate:
    {naive_ra_per100:.{d}e}
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
    samples = mgs_data.sample_attributes(
        bioproject, enrichment=Enrichment.VIRAL
    )

    for pathogen_name in ["sars_cov_2", "norovirus"]:
        pathogen = pathogens[pathogen_name]
        assert pathogen.pathogen_chars.na_type == NAType.RNA
        taxids = pathogen.pathogen_chars.taxids

        data = {
            "total_reads": np.array(
                [mgs_data.total_reads(bioproject)[s] for s in samples]
            ),
            "viral_reads": np.array(
                [mgs_data.viral_reads(bioproject, taxids)[s] for s in samples]
            ),
            "prevalence_per100k": np.array(
                stats.lookup_prevalence(samples, pathogen)
            ),
            "county": [s.county for s in samples.values()],
            "date": [s.date for s in samples.values()],
            "plant": [s.fine_location for s in samples.values()],
            "observation_type": "data",
        }

        fit = stats.fit_model(
            data=data,
            random_seed=1,
        )
        df = stats.fit_to_dataframe(fit, samples)
        df = pd.concat([pd.DataFrame(data), df], ignore_index=True)

        df.to_csv(
            f"fits/rothman-{pathogen_name}.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )

        naive_ra_per100 = (
            stats.per100k_to_per100
            * stats.naive_relative_abundance(
                data["viral_reads"],
                data["total_reads"],
                np.mean(data["prevalence_per100k"]),
            )
        )
        model_ra_per100 = pd.pivot_table(
            df, index="draws", values=["ra_per_one_percent"]
        )
        print_summary(pathogen_name, naive_ra_per100, model_ra_per100)


if __name__ == "__main__":
    start()
