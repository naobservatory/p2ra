from dataclasses import dataclass
from datetime import date
from typing import Any, TypeVar

import numpy as np
import pandas as pd
import stan  # type: ignore

from mgs import BioProject, Enrichment, MGSData, Sample, SampleAttributes
from pathogen_properties import Variable

per100k_to_per100 = 1e5 / 1e2


def is_match(
    sample_attrs: SampleAttributes,
    variable: Variable,
) -> bool:
    country, state, county = variable.get_location()
    start, end = variable.get_dates()
    assert isinstance(sample_attrs.date, date)
    return (
        (variable.taxid is None)  # TODO: allow other taxids
        and (country == sample_attrs.country)
        and ((state is None) or (state == sample_attrs.state))
        and ((county is None) or (county == sample_attrs.county))
        and (start <= sample_attrs.date <= end)
    )


V = TypeVar("V", bound=Variable)


def match_variables(
    samples: dict[Sample, SampleAttributes],
    vars: list[V],
) -> list[V]:
    keep = []
    for _, attrs in samples.items():
        matches = [v for v in vars if is_match(attrs, v)]
        # TODO: handle multiple matches
        assert len(matches) == 1
        keep.append(matches[0])
    return keep


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
    df["incidence_per100k"] = np.exp(df["theta"])
    df["ra_per_one_percent"] = per100k_to_per100 * np.exp(df["b"])
    df["observation_type"] = "posterior"
    return df


stan_code = """
data {
  int<lower=1> J;
  array[J] int<lower=0> y;
  vector[J] n;
  vector[J] x;
}
transformed data {
  vector[J] mu = log(x);
  real<lower=0> sigma = 0.5;
}
parameters {
  real b;                 // log conversion factor
  real<lower=0> phi;      // inverse overdispersion
  vector[J] theta;        // true log incidence
}
model {
  b ~ normal(0, 10);
  phi ~ gamma(2, 2);

  theta ~ normal(mu, sigma);
  y ~ neg_binomial_2_log(b + theta + log(n), phi);
}
generated quantities {
  array[J] int<lower=0> y_tilde
    = neg_binomial_2_log_rng(b + theta + log(n), phi);
}
"""


@dataclass
class Model:
    x: str
    data: dict[str, Any]
    samples: dict[Sample, SampleAttributes]
    y: str = "viral_reads"
    fit: None | stan.model.Model = None

    def fit_model(
        self, random_seed: int, num_chains: int = 4, num_samples: int = 1000
    ) -> None:
        stan_data = {
            "J": len(self.samples),
            "n": self.data["total_reads"],
            "x": self.data[self.x],
            "y": self.data[self.y],
        }
        model = stan.build(stan_code, data=stan_data, random_seed=random_seed)
        self.fit = model.sample(num_chains=num_chains, num_samples=num_samples)

    def get_fit_dataframe(self) -> pd.DataFrame:
        if not self.fit:
            raise ValueError("Model not fit yet")
        df = fit_to_dataframe(self.fit, self.samples)
        df = pd.concat([pd.DataFrame(self.data), df], ignore_index=True)
        return df


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
    incidences = pathogen.estimate_incidences()
    data = {
        "total_reads": np.array(
            [mgs_data.total_reads(bioproject)[s] for s in samples]
        ),
        "viral_reads": np.array(
            [mgs_data.viral_reads(bioproject, taxids)[s] for s in samples]
        ),
        "incidence_per100k": np.array(
            [
                inc.annual_infections_per_100k
                for inc in match_variables(samples, incidences)
            ]
        ),
        "county": [s.county for s in samples.values()],
        "date": [s.date for s in samples.values()],
        "plant": [s.fine_location for s in samples.values()],
        "observation_type": "data",
    }
    model = Model(data=data, samples=samples, x="incidence_per100k")
    model.fit_model(random_seed=random_seed)
    return model.get_fit_dataframe()
