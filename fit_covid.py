#!/usr/bin/env python3
from datetime import date

import numpy as np
import pandas as pd

import stats
from fit_rothman import per100k_to_per100, print_summary
from mgs import BioProject, Enrichment, MGSData, Sample, SampleAttributes
from pathogens import pathogens


def prevalence_by_state_county_date(
    pathogen: str,
) -> dict[tuple[str, str, date], float]:
    prevs = {}
    for estimate in pathogens[pathogen].estimate_prevalences():
        assert estimate.parsed_start == estimate.parsed_end
        country, state, county = estimate.target_location()
        assert country == "United States"
        assert state is not None
        assert county is not None
        est_date = estimate.parsed_start
        key = (state, county, est_date)
        assert key not in prevs
        prevs[key] = estimate.infections_per_100k
    return prevs


def lookup_prevalence(
    samples: dict[Sample, SampleAttributes], pathogen: str
) -> list[float]:
    lookup = prevalence_by_state_county_date(pathogen)
    prevs = []
    for _, attrs in samples.items():
        assert attrs.fine_location is not None
        assert attrs.county in [
            "Los Angeles County",
            "San Diego County",
            "Orange County",
        ]
        assert isinstance(attrs.date, date)
        prevs.append(lookup[("California", attrs.county, attrs.date)])
    return prevs


def fit_to_dataframe(
    fit, samples: dict[Sample, SampleAttributes]
) -> pd.DataFrame:
    df = pd.wide_to_long(
        fit.to_frame().reset_index(),
        stubnames=["y_tilde", "theta"],
        i="draws",
        j="sample",
        sep=".",
    ).reset_index()

    attrs = list(samples.values())

    def get_sample_attrs(attr: str):
        f = lambda i: getattr(attrs[i - 1], attr)
        return np.vectorize(f)

    df["date"] = get_sample_attrs("date")(df["sample"])
    df["county"] = get_sample_attrs("county")(df["sample"])
    df["plant"] = get_sample_attrs("fine_location")(df["sample"])
    df["total_reads"] = get_sample_attrs("reads")(df["sample"])

    df["viral_reads"] = df["y_tilde"]
    df["prevalence_per100k"] = np.exp(df["theta"])
    df["ra_per_one_percent"] = per100k_to_per100 * np.exp(df["b"])
    df["observation_type"] = "posterior"
    return df


def start():
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

    prevalence_per100k = np.array(lookup_prevalence(samples, pathogen))

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
    df = fit_to_dataframe(fit, samples)

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
    df = pd.concat([df, df_obs], ignore_index=True)

    df.to_csv(
        "fits/rothman-sars_cov_2.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    # TODO: Find a better way to get the once-per-draw stats
    model_ra_per100 = df[df["sample"] == 1]["ra_per_one_percent"]
    print_summary(pathogen, naive_ra_per100, model_ra_per100)


if __name__ == "__main__":
    start()
