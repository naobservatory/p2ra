from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Generic, TypeVar, Optional

import matplotlib  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
import numpy as np
import pandas as pd
import seaborn as sns  # type: ignore
import stan  # type: ignore
from scipy.stats import gamma, norm  # type: ignore

from mgs import BioProject, Enrichment, MGSData, Sample, SampleAttributes
from pathogen_properties import Predictor, Variable
from pathogens import pathogens

county_neighbors = {
    "Los Angeles County": [
        "Orange County",
        "San Diego County",
    ],
    "San Francisco County": [
        "Alameda County",
        "Marin County",
    ],
}

def county_is_close(county_a, county_b):
    def helper(county_1, county_2):
        return county_1 in county_neighbors and county_2 in county_neighbors[county_1]
    return helper(county_a, county_b) or helper(county_b, county_a)
            
def date_distance(start, end, target):
    if start <= target <= end:
        return 0
    return min(abs((start - target).days),
               abs((end - target).days))

def match_quality(
    sample_attrs: SampleAttributes,
    variable: Variable,
) -> Optional[int]:
    country, state, county = variable.get_location()
    start, end = variable.get_dates()
    assert isinstance(sample_attrs.date, date)

    if country != sample_attrs.country:
        return None
    if state is not None and state != sample_attrs.state:
        return None

    quality = 0
    if county is not None:
        if county == sample_attrs.county:
            # Prefer an exact match
            quality += 10
        elif not county_is_close(county, sample_attrs.county):
            # If no an exact match, require same metro area
            return None

    days_off = date_distance(start, end, sample_attrs.date)
    max_days_off = 7*2 # don't allow date matches more than two weeks out
    if days_off > max_days_off:
        return None
    quality -= days_off

    return quality

V = TypeVar("V", bound=Variable)


def lookup_variables(
    attrs: SampleAttributes,
    vars: list[V],
) -> list[V]:
    # Rank all matches by how close they are, then return all the ones tied for
    # best if there are any acceptable ones.
    #
    # We prefer matches that are temporally and geographically close.

    qualities = [
        (quality, var)
        for var in vars
        if (quality := match_quality(attrs, var)) is not None
    ]
    
    if not qualities:
        return []
    best_quality = max(quality for (quality, var) in qualities)
    return [var for (quality, var) in qualities if quality == best_quality]


P = TypeVar("P", bound=Predictor)


@dataclass
class DataPoint(Generic[P]):
    sample: Sample
    attrs: SampleAttributes
    viral_reads: int
    predictor: P


STANFILE = Path("model.stan")


@dataclass
class Model(Generic[P]):
    data: list[DataPoint[P]]
    fit: None | stan.fit.Fit = None
    dataframe: None | pd.DataFrame = None

    def fit_model(
        self, random_seed: int, num_chains: int = 4, num_samples: int = 1000
    ) -> None:
        with open(STANFILE, "r") as stanfile:
            stan_code = stanfile.read()
        stan_data = {
            "J": len(self.data),
            "y": np.array([dp.viral_reads for dp in self.data]),
            "n": np.array([dp.attrs.reads for dp in self.data]),
            "x": np.array([dp.predictor.get_data() for dp in self.data]),
        }
        model = stan.build(stan_code, data=stan_data, random_seed=random_seed)
        self.fit = model.sample(num_chains=num_chains, num_samples=num_samples)
        self.dataframe = self._make_dataframe()

    def _make_dataframe(self) -> pd.DataFrame:
        if self.fit is None:
            raise ValueError("Model not fit yet")

        df_input = pd.DataFrame(
            {
                "sample": [i + 1 for i, _ in enumerate(self.data)],
                "viral_reads": [dp.viral_reads for dp in self.data],
                "predictor": [dp.predictor.get_data() for dp in self.data],
                "observation_type": "data",
            }
        )

        df_output = pd.wide_to_long(
            self.fit.to_frame().reset_index(),
            stubnames=["y_tilde", "theta", "theta_std"],
            i="draws",
            j="sample",
            sep=".",
        ).reset_index()
        df_output["predictor"] = np.exp(df_output["theta"])
        df_output["ra_per_predictor"] = np.exp(df_output["b"])
        df_output["observation_type"] = "posterior"
        df_output.rename(columns={"y_tilde": "viral_reads"}, inplace=True)

        df = pd.concat([df_input, df_output], ignore_index=True)

        def get_sample_attrs(attr: str):
            f = lambda i: getattr(self.data[i - 1].attrs, attr)
            return np.vectorize(f)

        for attr in ["date", "county", "fine_location", "reads"]:
            df[attr] = get_sample_attrs(attr)(df["sample"])
        df.rename(columns={"reads": "total_reads"}, inplace=True)

        return df

    def get_per_draw_statistics(self) -> pd.DataFrame:
        if self.fit is None:
            raise ValueError("Model not fit yet")
        return self.fit.to_frame()

    def plot_posterior_histograms(self) -> matplotlib.figure.Figure:
        # TODO: Make sure this stays in sync with model.stan
        params = [
            ("phi", np.linspace(0, 6, 1000), gamma(2.0, scale=2.0)),
            ("b_std", np.linspace(-4, 4, 1000), norm(scale=2)),
        ]
        per_draw_df = self.get_per_draw_statistics()
        fig, axes = plt.subplots(
            1, len(params), layout="constrained", figsize=(6, 3)
        )
        for (param, x, prior), ax in zip(params, axes):
            posterior_hist(
                data=per_draw_df, param=param, prior_x=x, prior=prior, ax=ax
            )
        return fig

    def plot_posterior_samples(
        self, x: str, y: str, **kwargs
    ) -> sns.FacetGrid:
        # Plot posterior predictive draws
        data = self.dataframe
        if data is None:
            raise ValueError("Model not fit yet")
        g = sns.relplot(
            data=data[
                (data["observation_type"] == "posterior") & (data["draws"] < 9)
            ],
            x=x,
            y=y,
            col="draws",
            col_wrap=3,
            height=4,
            **kwargs,
        )
        g.set_titles("Posterior draw {col_name:1.0f}")
        # Plot data
        ax = g.facet_axis(0, 0)
        for col in ax.collections:
            col.remove()
        sns.scatterplot(
            data=data[data["observation_type"] == "data"],
            x=x,
            y=y,
            ax=ax,
            legend=False,
            **kwargs,
        )
        ax.set_title("Observed", fontdict={"size": 10})
        return g


def choose_predictor(predictors: list[Predictor]) -> Predictor:
    # TODO: allow other taxids
    non_taxid_predictors = [
        predictor for predictor in predictors if not predictor.taxid]
    assert len(non_taxid_predictors) == 1
    return non_taxid_predictors[0]


def build_model(
    mgs_data: MGSData,
    bioproject: BioProject,
    pathogen_name: str,
    predictor: str,
) -> Model:
    pathogen = pathogens[pathogen_name]
    taxids = pathogen.pathogen_chars.taxids
    samples = mgs_data.sample_attributes(
        bioproject, enrichment=Enrichment.VIRAL
    )
    predictors: list[Predictor]
    if predictor == "incidence":
        predictors = pathogen.estimate_incidences()
    elif predictor == "prevalence":
        predictors = pathogen.estimate_prevalences()
    else:
        raise ValueError(
            f"{predictor} must be one of 'incidence' or 'prevalence'"
        )
    data = [
        DataPoint(
            sample=s,
            attrs=attrs,
            viral_reads=mgs_data.viral_reads(bioproject, taxids)[s],
            predictor=choose_predictor(lookup_variables(attrs, predictors)),
        )
        for s, attrs in samples.items()
    ]
    return Model(data=data)


def posterior_hist(data, param: str, prior_x, prior, ax=None):
    sns.lineplot(
        x=prior_x, y=prior.pdf(prior_x), color="black", label="prior", ax=ax
    )
    sns.histplot(data=data, x=param, stat="density", bins=40, ax=ax)
    return ax
