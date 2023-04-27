#!/usr/bin/env python3

import numpy as np  # type: ignore

import stats
from mgs import BioProject, MGSData
from pathogens import pathogens

if __name__ == "__main__":
    bioproject = BioProject("PRJNA729801")  # Rothman
    pathogen = "norovirus"
    taxid = pathogens[pathogen].pathogen_chars.taxid

    mgs_data = MGSData.from_repo()
    samples = mgs_data.sample_attributes(bioproject).keys()
    all_reads = [mgs_data.total_reads(bioproject)[s] for s in samples]
    virus_reads = [mgs_data.viral_reads(bioproject, taxid)[s] for s in samples]

    prevalence_estimates = pathogens[pathogen].estimate_prevalences()
    assert len(prevalence_estimates) == 1
    prevalence_per100k = prevalence_estimates[0].infections_per_100k

    per100k_to_per100 = 1e5 / 1e2

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
    model_ra_per100 = per100k_to_per100 * np.exp(fit["b"])
    percentiles = [5, 25, 50, 75, 95]
    sep = "\t\t"
    title = f"Relative abundance of {pathogen} at 1% prevalence"
    print("-" * len(title))
    print(title)
    print("-" * len(title))
    print("Naive estimate:")
    print(f"{naive_ra_per100:0.5f}")
    print("Posterior arithmetic mean:")
    print(f"{np.mean(model_ra_per100):0.5f}")
    print("Posterior geometric mean:")
    print(f"{np.exp(np.mean(np.log(model_ra_per100))):0.5f}")
    print("Posterior quantiles:")
    print(*(f"{p}%" for p in percentiles), sep=sep)
    print(
        *(f"{x:0.5f}" for x in np.percentile(model_ra_per100, percentiles)),
        sep=sep,
    )
