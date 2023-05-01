#!/usr/bin/env python3
from textwrap import dedent

import numpy as np

import stats
from mgs import BioProject, MGSData
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


per100k_to_per100 = 1e5 / 1e2

if __name__ == "__main__":
    bioproject = BioProject("PRJNA729801")  # Rothman

    mgs_data = MGSData.from_repo()
    samples = mgs_data.sample_attributes(bioproject).keys()
    all_reads = [mgs_data.total_reads(bioproject)[s] for s in samples]

    for pathogen_name, pathogen in pathogens.items():
        taxid = pathogen.pathogen_chars.taxid
        virus_reads = [
            mgs_data.viral_reads(bioproject, taxid)[s] for s in samples
        ]

        prevalence_estimates = pathogen.estimate_prevalences()
        prevalence_per100k = np.mean(
            [est.infections_per_100k for est in prevalence_estimates]
        )

        naive_ra_per100 = per100k_to_per100 * stats.naive_relative_abundance(
            virus_reads,
            all_reads,
            prevalence_per100k,
        )

        fit = stats.fit_model(
            num_samples=len(samples),
            viral_read_counts=virus_reads,
            total_read_counts=all_reads,
            mean_log_prevalence=np.log(prevalence_per100k),
            std_log_prevalence=1.0,
            random_seed=1,
        )
        # TODO: Wrap the model fit so that we aren't exposed to stan variables
        model_ra_per100 = per100k_to_per100 * np.exp(fit["b"])
        print_summary(pathogen_name, naive_ra_per100, model_ra_per100)
