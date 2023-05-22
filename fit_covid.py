#!/usr/bin/env python3
from datetime import date

import numpy as np
import pandas as pd

import stats
from fit_rothman import per100k_to_per100, print_summary
from mgs import BioProject, Enrichment, MGSData, Sample, SampleAttributes
from pathogen_properties import NAType, Prevalence
from pathogens import pathogens


def is_match(
    prevalence: Prevalence,
    sample_attrs: SampleAttributes,
) -> bool:
    country, state, county = prevalence.get_location()
    start, end = prevalence.get_dates()
    assert isinstance(sample_attrs.date, date)
    return (
        (prevalence.taxid is None)  # TODO: allow other taxids
        and (country == sample_attrs.country)
        and ((state is None) or (state == sample_attrs.state))
        and ((county is None) or (county == sample_attrs.county))
        and (start <= sample_attrs.date <= end)
    )


def lookup_prevalence(
    samples: dict[Sample, SampleAttributes],
    pathogen,
) -> list[float]:
    prev_estimates = pathogen.estimate_prevalences()
    prevs = []
    for _, attrs in samples.items():
        matches = [p for p in prev_estimates if is_match(p, attrs)]
        # TODO: handle multiple matches
        assert len(matches) == 1
        prevs.append(matches[0].infections_per_100k)
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

    for pathogen_name in ["sars_cov_2", "norovirus"]:
        pathogen = pathogens[pathogen_name]
        taxids = pathogen.pathogen_chars.taxids

        data = {
            "total_reads": np.array(
                [mgs_data.total_reads(bioproject)[s] for s in samples]
            ),
            "viral_reads": np.array(
                [mgs_data.viral_reads(bioproject, taxids)[s] for s in samples]
            ),
            "prevalence_per100k": np.array(
                lookup_prevalence(samples, pathogen)
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
        df = fit_to_dataframe(fit, samples)
        df = pd.concat([pd.DataFrame(data), df], ignore_index=True)

        df.to_csv(
            f"fits/rothman-{pathogen_name}.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )

        naive_ra_per100 = per100k_to_per100 * stats.naive_relative_abundance(
            data["viral_reads"],
            data["total_reads"],
            np.mean(data["prevalence_per100k"]),
        )
        model_ra_per100 = pd.pivot_table(
            df, index="draws", values=["ra_per_one_percent"]
        )
        print_summary(pathogen_name, naive_ra_per100, model_ra_per100)


if __name__ == "__main__":
    start()
