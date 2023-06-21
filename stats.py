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
    predictor: P


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
    fine_locations: list[str | None] = field(init=False)
    fit: None | stan.fit.Fit = None
    dataframe: None | pd.DataFrame = None

    def __post_init__(self) -> None:
        with open(STANFILE, "r") as stanfile:
            stan_code = Template(stanfile.read()).substitute(**HYPERPARAMS)
        print(stan_code)
        # TODO: Make it more automatic to associate fine locations with coeffs
        self.fine_locations = sorted(
            list(set(dp.attrs.fine_location for dp in self.data)), key=str
        )
        stan_data = {
            "J": len(self.data),
            "y": np.array([dp.viral_reads for dp in self.data]),
            "n": np.array([dp.attrs.reads for dp in self.data]),
            "x": np.array([dp.predictor.get_data() for dp in self.data]),
            "L": len(self.fine_locations),
            "ll": [
                # Stan vectors are one-indexed
                self.fine_locations.index(dp.attrs.fine_location) + 1
                for dp in self.data
            ],
        }
        self.model = stan.build(
            stan_code, data=stan_data, random_seed=self.random_seed
        )

    def fit_model(self, num_chains: int = 4, num_samples: int = 1000) -> None:
        self.fit = self.model.sample(
            num_chains=num_chains, num_samples=num_samples
        )
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

    def plot_data_scatter(self) -> matplotlib.figure.Figure:
        df = pd.DataFrame(
            {
                "county": [dp.attrs.county for dp in self.data],
                "fine_location": [dp.attrs.fine_location for dp in self.data],
                "relative_abundance": np.array(
                    [dp.viral_reads for dp in self.data]
                )
                / np.array([dp.attrs.reads for dp in self.data]),
                "predictor": np.array(
                    [dp.predictor.get_data() for dp in self.data]
                ),
            }
        )
        fig, ax = plt.subplots(1, 1)
        sns.scatterplot(
            data=df,
            x="predictor",
            y="relative_abundance",
            ax=ax,
            style="county",
            hue="fine_location",
            hue_order=self.fine_locations,
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
        per_draw_df = self.get_per_draw_statistics()
        fig, axes = plt.subplots(
            1, len(params), layout="constrained", figsize=(6, 3)
        )
        for (param, x, prior), ax in zip(params, axes):
            posterior_hist(
                data=per_draw_df, param=param, prior_x=x, prior=prior, ax=ax
            )
        return fig

    def plot_violin(self) -> matplotlib.figure.Figure:
        per_draw_df = self.get_per_draw_statistics()
        df = pd.DataFrame()
        # TODO: this probably goes elsewhere
        for i, loc in enumerate(self.fine_locations):
            df[loc] = per_draw_df[f"b_l.{i+1}"]
        df["Overall"] = per_draw_df["mu"]
        fig, ax = plt.subplots(1, 1)
        sns.violinplot(data=df, ax=ax)
        ax.set_ylabel("Standardized coefficient")
        ax.set_xlabel("Sampling location")
        return fig

    def plot_joint_posterior(self, x: str, y: str) -> matplotlib.figure.Figure:
        per_draw_df = self.get_per_draw_statistics()
        fig, ax = plt.subplots(1, 1)
        sns.kdeplot(
            data=per_draw_df,
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

    def plot_figures(self, path: Path, prefix: str) -> None:
        assert self.fit is not None
        data_scatter = self.plot_data_scatter()
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
                style="county",
                hue="fine_location",
                hue_order=self.fine_locations,
            )
            if y == "predictor":
                g.set(yscale="log")
            g.savefig(path / f"{prefix}-{y}_vs_{x}.pdf")
        plt.close("all")


def choose_predictor(predictors: list[Predictor]) -> Predictor:
    assert len(predictors) == 1
    return predictors[0]


def build_model(
    mgs_data: MGSData,
    bioproject: BioProject,
    predictors: list[Predictor],
    taxids: frozenset[TaxID],
    random_seed: int,
) -> Model:
    samples = mgs_data.sample_attributes(
        bioproject, enrichment=Enrichment.VIRAL
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
    return Model(data=data, random_seed=random_seed)


def posterior_hist(data, param: str, prior_x, prior, ax=None):
    sns.lineplot(
        x=prior_x, y=prior.pdf(prior_x), color="black", label="prior", ax=ax
    )
    sns.histplot(data=data, x=param, stat="density", bins=40, ax=ax)
    return ax
