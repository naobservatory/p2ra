#!/usr/bin/env python3
from textwrap import dedent

import numpy as np
import pandas as pd

import stats
from fit_rothman import per100k_to_per100, print_summary
from mgs import BioProject, Enrichment, MGSData
from pathogens import pathogens

if __name__ == "__main__":
    bioproject = BioProject("PRJNA729801")  # Rothman

    mgs_data = MGSData.from_repo()
    samples = mgs_data.sample_attributes(
        bioproject, enrichment=Enrichment.VIRAL
    )
    all_reads = np.array(
        [mgs_data.total_reads(bioproject)[s] for s in samples]
    )

    pathogen = "sars_cov_2"
    taxids = pathogens[pathogen].pathogen_chars.taxids
    virus_reads = np.array(
        [mgs_data.viral_reads(bioproject, taxids)[s] for s in samples]
    )

    prevalence_by_loc_date = {}
    for estimate in pathogens[pathogen].estimate_prevalences():
        assert estimate.parsed_start == estimate.parsed_end
        key = (estimate.county, estimate.parsed_start)
        assert key not in prevalence_by_loc_date
        prevalence_by_loc_date[key] = estimate.infections_per_100k

    prevalence_per100k = np.zeros(len(samples))
    for i, (sample, attrs) in enumerate(samples.items()):
        assert attrs.fine_location is not None
        assert attrs.county is not None
        prevalence_per100k[i] = prevalence_by_loc_date[
            (attrs.county, attrs.date)
        ]

    naive_ra_per100 = per100k_to_per100 * stats.naive_relative_abundance(
        virus_reads,
        all_reads,
        np.mean(prevalence_per100k),
    )

    fit = stats.fit_model(
        num_samples=len(virus_reads),
        viral_read_counts=virus_reads,
        total_read_counts=all_reads,
        mean_log_prevalence=np.log(prevalence_per100k),
        std_log_prevalence=0.5,
        random_seed=1,
    )

    # TODO: do this more neatly
    df_obs = pd.DataFrame(
        {
            "viral_reads": virus_reads,
            "total_reads": all_reads,
            "prevalence_per100k": prevalence_per100k,
            "county": [s.county for s in samples.values()],
            "date": [s.date for s in samples.values()],
            "plant": [s.fine_location for s in samples.values()],
            "observation_type": "data",
        }
    )
    df_obs.to_csv("covid_input.tsv", sep="\t")

    # TODO: Wrap the model fit so that we aren't exposed to stan variables
    model_ra_per100 = per100k_to_per100 * np.exp(fit["b"])
    print_summary(pathogen, naive_ra_per100, model_ra_per100)
    pp_virus_reads = fit["y_tilde"]
    print(np.mean(virus_reads))
    means = np.mean(pp_virus_reads, axis=0)
    d = 1
    sep = " " * 4
    percentiles = [5, 25, 50, 75, 95]
    percentile_values = np.percentile(means, percentiles)
    print(
        dedent(
            f"""Posterior quantiles:
    {sep.join(f"{p:>{d+5}}%" for p in percentiles)}
    {sep.join(f"{x:.{d}e}" for x in percentile_values)}"""
        )
    )
    print(np.max(virus_reads))
    maxes = np.max(pp_virus_reads, axis=0)
    percentile_values = np.percentile(maxes, percentiles)
    print(
        dedent(
            f"""Posterior quantiles:
    {sep.join(f"{p:>{d+5}}%" for p in percentiles)}
    {sep.join(f"{x:.{d}e}" for x in percentile_values)}"""
        )
    )
    stats.save_fit(fit, "fits/rothman-sars_cov_2.tsv")
