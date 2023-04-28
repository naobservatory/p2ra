#!/usr/bin/env python3
import numpy as np

import stats
from fit_rothman import per100k_to_per100, print_summary
from mgs import BioProject, MGSData
from pathogens import pathogens

plant_counties = {
    "HTP": "Los Angeles",
    "SJ": "Los Angeles",
    "JWPCP": "Los Angeles",
    "OC": "Orange",
    "PL": "San Diego",
    "SB": "San Diego",
    "NC": "San Diego",
}

if __name__ == "__main__":
    bioproject = BioProject("PRJNA729801")  # Rothman

    mgs_data = MGSData.from_repo()
    samples = mgs_data.sample_attributes(bioproject)
    all_reads = np.array(
        [mgs_data.total_reads(bioproject)[s] for s in samples]
    )

    pathogen = "sars_cov_2"
    taxid = pathogens[pathogen].pathogen_chars.taxid
    virus_reads = np.array(
        [mgs_data.viral_reads(bioproject, taxid)[s] for s in samples]
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
        try:
            county = plant_counties[attrs.fine_location]
        except KeyError:
            prevalence_per100k[i] = np.nan
        else:
            prevalence_per100k[i] = prevalence_by_loc_date[
                (county, attrs.date)
            ]

    virus_reads = virus_reads[~np.isnan(prevalence_per100k)]
    all_reads = all_reads[~np.isnan(prevalence_per100k)]
    prevalence_per100k = prevalence_per100k[~np.isnan(prevalence_per100k)]

    naive_ra_per100 = per100k_to_per100 * stats.naive_relative_abundance(
        virus_reads,
        all_reads,
        # Stop mypy from being persnickety while I sort out numpy annotations
        float(np.mean(prevalence_per100k)),
    )

    fit = stats.fit_model(
        num_samples=len(virus_reads),
        viral_read_counts=virus_reads,
        total_read_counts=all_reads,
        mean_log_prevalence=np.log(prevalence_per100k),
        std_log_prevalence=0.5,
        random_seed=1,
    )
    # TODO: Wrap the model fit so that we aren't exposed to stan variables
    model_ra_per100 = per100k_to_per100 * np.exp(fit["b"])
    print_summary(pathogen, naive_ra_per100, model_ra_per100)
