from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from string import Template
from typing import Generic, Optional, TypeVar

import matplotlib  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
import numpy as np
import pandas as pd
import seaborn as sns  # type: ignore
import stan  # type: ignore
from scipy.stats import gamma, norm  # type: ignore

from mgs import BioProject, Enrichment, MGSData, Sample, SampleAttributes
from pathogen_properties import Predictor, TaxID, Variable

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
        return (
            county_1 in county_neighbors
            and county_2 in county_neighbors[county_1]
        )

    return helper(county_a, county_b) or helper(county_b, county_a)


def date_distance(start, end, target):
    if start <= target <= end:
        return 0
    return min(abs((start - target).days), abs((end - target).days))


def match_quality(
    sample_attrs: SampleAttributes,
    variable: Variable,
) -> Optional[int]:
    country, state, county = variable.get_location()
    start, end = variable.get_dates()
    assert isinstance(sample_attrs.date, date)

    if country != sample_attrs.country:
        return None

    quality = 0
    if state is not None:
        if state != sample_attrs.state:
            return None
        # Prefer the specific state.
        quality += 20

    if county is not None:
        if county == sample_attrs.county:
            # Prefer an exact match
            quality += 10
        elif not county_is_close(county, sample_attrs.county):
            # If no an exact match, require same metro area
            return None

    days_off = date_distance(start, end, sample_attrs.date)
    max_days_off = 7 * 2  # don't allow date matches more than two weeks out
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
    best_quality = max(quality for (quality, _) in qualities)
    return [var for (quality, var) in qualities if quality == best_quality]


P = TypeVar("P", bound=Predictor)


@dataclass
class DataPoint(Generic[P]):
    sample: Sample
    attrs: SampleAttributes
    viral_reads: int
    predictor: P | None

    def get_predictor_value(self) -> float:
        # TODO: If we update the Stan code to allow some samples to be
        # missing predictors, have this return the value Stan expects for
        # missing input data (e.g. NaN)
        if self.predictor is None:
            raise NotImplementedError(
                f"Data point for sample {self.sample} missing predictor"
            )
        else:
            return self.predictor.get_data()


STANFILE = Path("model.stan")

# TODO: Make this configurable
HYPERPARAMS = {
    "mu_sigma": 4,
    "sigma_alpha": 2,
    "sigma_beta": 1,
    "tau_alpha": 2,
    "tau_beta": 1,
}


@dataclass
class Model(Generic[P]):
    data: list[DataPoint[P]]
    random_seed: int
    model: stan.model.Model = field(init=False)
    locations: list[str | None] = field(init=False)
    input_df: pd.DataFrame = field(init=False)
    fit: None | stan.fit.Fit = None
    output_df: None | pd.DataFrame = None

    def __post_init__(self) -> None:
        with open(STANFILE, "r") as stanfile:
            stan_code = Template(stanfile.read()).substitute(**HYPERPARAMS)
        self.input_df = pd.DataFrame(
            {
                # Stan vectors are 1-indexed
                "sample": [i + 1 for i, _ in enumerate(self.data)],
                "viral_reads": [dp.viral_reads for dp in self.data],
                "total_reads": [dp.attrs.reads for dp in self.data],
                "predictor": [dp.get_predictor_value() for dp in self.data],
                "fine_location": [dp.attrs.fine_location for dp in self.data],
                "date": [dp.attrs.date for dp in self.data],
                "county": [dp.attrs.county for dp in self.data],
                "relative_abundance": np.array(
                    [dp.viral_reads for dp in self.data]
                )
                / np.array([dp.attrs.reads for dp in self.data]),
            }
        )
        # TODO: Make it more automatic to associate fine locations with coeffs
        self.locations = sorted(
            list(set(dp.attrs.fine_location for dp in self.data)), key=str
        ) + ["Overall"]
        stan_data = {
            "J": len(self.data),
            "y": self.input_df.viral_reads.to_numpy(),
            "n": self.input_df.total_reads.to_numpy(),
            "x": self.input_df.predictor.to_numpy(),
            # Overall is not a location
            "L": len(self.locations) - 1,
            "ll": [
                # Stan vectors are one-indexed
                self.locations.index(loc) + 1
                for loc in self.input_df.fine_location
            ],
        }
        self.model = stan.build(
            stan_code, data=stan_data, random_seed=self.random_seed
        )

    def fit_model(self, num_chains: int = 4, num_samples: int = 1000) -> None:
        self.fit = self.model.sample(
            num_chains=num_chains, num_samples=num_samples
        )
        self.output_df = self.fit.to_frame()

    def get_output_by_sample(self) -> pd.DataFrame:
        if self.output_df is None:
            raise ValueError("Model not fit yet")

        df = pd.wide_to_long(
            self.output_df.reset_index(),
            stubnames=["y_tilde", "theta", "theta_std"],
            i="draws",
            j="sample",
            sep=".",
        ).reset_index()
        df["predictor"] = np.exp(df["theta"])
        df.rename(columns={"y_tilde": "viral_reads"}, inplace=True)

        def get_sample_attrs(attr: str):
            f = lambda i: getattr(self.data[i - 1].attrs, attr)
            return np.vectorize(f)

        for attr in ["date", "county", "fine_location", "reads"]:
            df[attr] = get_sample_attrs(attr)(df["sample"])
        df.rename(columns={"reads": "total_reads"}, inplace=True)

        return df

    def get_coefficients(self) -> pd.DataFrame:
        if self.output_df is None:
            raise ValueError("Model not fit yet")
        cols = ["b", "ra_at_1in100"]
        coeffs = pd.wide_to_long(
            self.output_df.reset_index(),
            stubnames=cols,
            i="draws",
            j="location_index",
            sep=".",
        ).reset_index()
        coeffs["location"] = np.array(self.locations)[
            coeffs["location_index"] - 1
        ]
        return coeffs[["location"] + cols]

    def plot_data_scatter(self, **kwargs) -> matplotlib.figure.Figure:
        fig, ax = plt.subplots(1, 1)
        sns.scatterplot(
            data=self.input_df,
            x="predictor",
            y="relative_abundance",
            ax=ax,
            hue="fine_location",
            hue_order=self.locations,
            **kwargs,
        )
        sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
        return fig

    def plot_posterior_histograms(self) -> matplotlib.figure.Figure:
        # TODO: Make sure this stays in sync with model.stan
        params = [
            (
                "sigma",
                np.linspace(0, 6, 1000),
                gamma(
                    HYPERPARAMS["sigma_alpha"],
                    scale=1 / HYPERPARAMS["sigma_beta"],
                ),
            ),
            (
                "mu",
                np.linspace(-8, 4, 1000),
                norm(scale=HYPERPARAMS["mu_sigma"]),
            ),
            (
                "tau",
                np.linspace(0, 6, 1000),
                gamma(
                    HYPERPARAMS["tau_alpha"],
                    scale=1 / HYPERPARAMS["tau_beta"],
                ),
            ),
        ]
        fig, axes = plt.subplots(
            1, len(params), layout="constrained", figsize=(6, 3)
        )
        for (param, x, prior), ax in zip(params, axes):
            posterior_hist(
                data=self.output_df, param=param, prior_x=x, prior=prior, ax=ax
            )
        return fig

    def plot_violin(self) -> matplotlib.figure.Figure:
        fig, ax = plt.subplots(1, 1)
        sns.violinplot(
            data=self.get_coefficients(), x="location", y="b", ax=ax
        )
        ax.set_ylabel("Standardized coefficient")
        ax.set_xlabel("Sampling location")
        return fig

    def plot_joint_posterior(self, x: str, y: str) -> matplotlib.figure.Figure:
        fig, ax = plt.subplots(1, 1)
        sns.kdeplot(
            data=self.output_df,
            ax=ax,
            x=x,
            y=y,
            fill=True,
            levels=100,
            cmap="mako",
            cbar=True,
        )
        return fig

    def plot_posterior_samples(
        self, x: str, y: str, **kwargs
    ) -> sns.FacetGrid:
        # Plot posterior predictive draws
        if self.output_df is None:
            raise ValueError("Model not fit yet")
        posterior_draws = self.get_output_by_sample()
        g = sns.relplot(
            data=posterior_draws[posterior_draws.draws < 9],
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
            data=self.input_df,
            x=x,
            y=y,
            ax=ax,
            legend=False,
            **kwargs,
        )
        ax.set_title("Observed", fontdict={"size": 10})
        return g

    def plot_figures(self, path: Path, prefix: str) -> None:
        assert self.fit is not None
        if any(self.input_df["county"]):
            style = "county"
        else:
            style = None
        data_scatter = self.plot_data_scatter(style=style)
        data_scatter.savefig(
            path / f"{prefix}-datascatter.pdf", bbox_inches="tight"
        )
        fig_hist = self.plot_posterior_histograms()
        fig_hist.savefig(path / f"{prefix}-posthist.pdf")
        fig_viol = self.plot_violin()
        fig_viol.savefig(path / f"{prefix}-violin.pdf")
        for x, y in [("mu", "sigma"), ("mu", "tau"), ("sigma", "tau")]:
            fig = self.plot_joint_posterior(x, y)
            fig.savefig(path / f"{prefix}-{y}_vs_{x}.pdf")
        for x, y in [
            ("date", "viral_reads"),
            ("date", "predictor"),
            ("predictor", "viral_reads"),
        ]:
            g = self.plot_posterior_samples(
                x,
                y,
                style=style,
                hue="fine_location",
                hue_order=self.locations,
            )
            if y == "predictor":
                g.set(yscale="log")
            g.savefig(path / f"{prefix}-{y}_vs_{x}.pdf")
        plt.close("all")


def choose_predictor(predictors: list[Predictor]) -> Predictor | None:
    if len(predictors) == 0:
        return None
    elif len(predictors) == 1:
        return predictors[0]
    else:
        raise NotImplementedError("More than one matching predictor")


def build_model(
    mgs_data: MGSData,
    bioprojects: list[BioProject],
    predictors: list[Predictor],
    taxids: frozenset[TaxID],
    random_seed: int,
    enrichment: Optional[Enrichment],
) -> Model | None:
    sample_attributes = {}  # sample -> attributes
    study_viral_reads = {}  # sample -> viral_reads
    for bioproject in bioprojects:
        sample_attributes.update(
            mgs_data.sample_attributes(bioproject, enrichment=enrichment)
        )
        study_viral_reads.update(mgs_data.viral_reads(bioproject, taxids))
    data = [
        DataPoint(
            sample=sample,
            attrs=attrs,
            viral_reads=study_viral_reads[sample],
            predictor=choose_predictor(lookup_variables(attrs, predictors)),
        )
        for sample, attrs in sample_attributes.items()
    ]
    # No predictors found
    if all(point.predictor is None for point in data):
        return None
    else:
        return Model(data=data, random_seed=random_seed)


def posterior_hist(data, param: str, prior_x, prior, ax=None):
    sns.lineplot(
        x=prior_x, y=prior.pdf(prior_x), color="black", label="prior", ax=ax
    )
    sns.histplot(data=data, x=param, stat="density", bins=40, ax=ax)
    return ax
