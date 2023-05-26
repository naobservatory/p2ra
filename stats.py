from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Generic, TypeVar

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


def lookup_variable(
    attrs: SampleAttributes,
    vars: list[V],
) -> V:
    matches = [v for v in vars if is_match(attrs, v)]
    assert len(matches) == 1
    return matches[0]


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

    def get_per_draw_statistics(
        self, values: list[str] = ["phi", "b", "b_std", "ra_per_predictor"]
    ) -> pd.DataFrame:
        if self.dataframe is None:
            raise ValueError("Model not fit yet")
        return pd.pivot_table(
            self.dataframe,
            index="draws",
            values=values,
        )

    def plot_posterior_histograms(self) -> matplotlib.figure.Figure:
        # TODO: Make sure this stays in sync with model.stan
        params = [
            ("phi", np.linspace(0, 6, 1000), gamma(2.0, scale=2.0)),
            ("b_std", np.linspace(-4, 4, 1000), norm(scale=2)),
        ]
        per_draw_df = self.get_per_draw_statistics([p for p, _, _ in params])
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
            predictor=lookup_variable(attrs, predictors),
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
