#!/usr/bin/env python3

from collections import Counter

import numpy as np  # type: ignore

from mgs import MGSData, count_reads
import stats
from pathogens import pathogens

if __name__ == "__main__":
    bioproject = "PRJNA729801"  # Rothman
    pathogen = "norovirus"

    mgs_data = MGSData.from_repo()
    samples = mgs_data.samples_by_bioproject(bioproject)

    sample_attribs = {s: mgs_data.samples[s] for s in samples}

    all_reads: Counter[str] = Counter()
    for attribs in sample_attribs.values():
        assert attribs.fine_location is not None
        all_reads[attribs.fine_location] += attribs.reads

    studies = all_reads.keys()
    samples_by_study = {
        study: [s for s in samples if sample_attribs[s].fine_location == study]
        for study in studies
    }

    taxid = pathogens[pathogen].pathogen_chars.taxid
    subtree = mgs_data.tax_tree[taxid]
    assert subtree is not None

    virus_counts = count_reads(subtree, mgs_data.read_counts)
    virus_reads = {
        study: sum(virus_counts[s] for s in samples_by_study[study])
        for study in studies
    }

    estimates = pathogens[pathogen].estimate_prevalences()
    prevalence_per100k = estimates[0].infections_per_100k

    naive_ra_per100 = (1e5 / 1e2) * stats.naive_ra(
        [virus_reads[s] for s in studies],
        [all_reads[s] for s in studies],
        prevalence_per100k,
    )

    fit = stats.fit_model(
        num_studies=len(studies),
        viral_read_counts=[virus_reads[s] for s in studies],
        total_read_counts=[all_reads[s] for s in studies],
        mean_log_prevalence=np.log(prevalence_per100k),
        std_log_prevalence=1.0,
        random_seed=1,
    )
    model_ra_per100 = np.exp(fit["b"]) * (1e5 / 1e2)
    percentiles = [5, 25, 50, 75, 95]
    sep = "\t\t"
    title = "Relative abundance of norovirus at 1% prevalence"
    print("-" * len(title))
    print(title)
    print("-" * len(title))
    print("Naive estimate:")
    print(f"{naive_ra_per100:0.5f}")
    print("Posterior geometric mean:")
    print(f"{np.exp(np.mean(np.log(model_ra_per100))):0.5f}")
    print("Posterior quantiles:")
    print(*(f"{p}%" for p in percentiles), sep=sep)
    print(
        *(f"{x:0.5f}" for x in np.percentile(model_ra_per100, percentiles)),
        sep=sep,
    )
