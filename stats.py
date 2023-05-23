from datetime import date

import numpy as np
import numpy.typing as npt
import pandas as pd
import stan  # type: ignore

from mgs import BioProject, Enrichment, MGSData, Sample, SampleAttributes
from pathogen_properties import Prevalence

per100k_to_per100 = 1e5 / 1e2


def naive_relative_abundance(
    virus_counts: npt.ArrayLike,
    all_counts: npt.ArrayLike,
    prev_per_100k: np.float_,
) -> float:
    total_virus = np.sum(virus_counts)
    total_counts = np.sum(all_counts)
    return total_virus / total_counts / prev_per_100k


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


stan_code = """
data {
  int<lower=1> J;
  array[J] int<lower=0> viral_reads;
  vector[J] total_reads;
  vector[J] prevalence_per100k;
}
transformed data {
  vector[J] mu = log(prevalence_per100k);
  real<lower=0> sigma = 0.5;
}
parameters {
  real b;                 // log conversion factor
  real<lower=0> phi;      // inverse overdispersion
  vector[J] theta;        // true log prevalence
}
model {
  b ~ normal(0, 10);
  phi ~ gamma(2, 2);

  theta ~ normal(mu, sigma);
  viral_reads ~ neg_binomial_2_log(b + theta + log(total_reads), phi);
}
generated quantities {
  array[J] int<lower=0> y_tilde
    = neg_binomial_2_log_rng(b + theta + log(total_reads), phi);
}
"""


def fit_model(
    mgs_data: MGSData,
    bioproject: BioProject,
    pathogen,
    random_seed: int,
) -> pd.DataFrame:
    samples = mgs_data.sample_attributes(
        bioproject, enrichment=Enrichment.VIRAL
    )
    taxids = pathogen.pathogen_chars.taxids
    data = {
        "total_reads": np.array(
            [mgs_data.total_reads(bioproject)[s] for s in samples]
        ),
        "viral_reads": np.array(
            [mgs_data.viral_reads(bioproject, taxids)[s] for s in samples]
        ),
        "prevalence_per100k": np.array(lookup_prevalence(samples, pathogen)),
        "county": [s.county for s in samples.values()],
        "date": [s.date for s in samples.values()],
        "plant": [s.fine_location for s in samples.values()],
        "observation_type": "data",
    }
    stan_data = {
        k: v
        for k, v in data.items()
        if k in ["viral_reads", "total_reads", "prevalence_per100k"]
    }
    stan_data["J"] = len(samples)
    model = stan.build(stan_code, data=stan_data, random_seed=random_seed)
    fit = model.sample(num_chains=4, num_samples=1000)
    df = fit_to_dataframe(fit, samples)
    df = pd.concat([pd.DataFrame(data), df], ignore_index=True)
    return df
