#!/usr/bin/env python3
from pathlib import Path

import matplotlib.pyplot as plt  # type: ignore
import pandas as pd
import seaborn as sns  # type: ignore

from mgs import rna_bioprojects

study_order = ["crits_christoph", "rothman", "spurbeck"]


def plot_overall(data: pd.DataFrame, input_data: pd.DataFrame) -> plt.Figure:
    viral_reads = count_viral_reads(input_data)
    plotting_order = viral_reads.sort_values(
        ["reads_by_pathogen", "pathogen"]
    ).reset_index()
    fig = plt.figure(figsize=(6, 10))
    ax = fig.add_subplot(1, 1, 1)
    sns.boxenplot(
        ax=ax,
        data=data[data.location == "Overall"],
        x="ra_at_1in1000",
        y="pathogen",
        order=plotting_order.pathogen.unique(),
        hue="study",
        hue_order=study_order,
        showfliers=False,
    )
    for num_reads, artist_list in zip(
        plotting_order.viral_reads, ax.collections
    ):
        artist_list.set_alpha(min(num_reads / 10 + 0.02, 1.0))
    ax.set_xscale("log")
    ax.set_ylabel("")
    return fig


def plot_by_location(
    data: pd.DataFrame, input_data: pd.DataFrame
) -> plt.Figure:
    viral_reads = count_viral_reads(input_data, by_location=True)
    plotting_order = viral_reads.sort_values(
        ["reads_by_pathogen", "pathogen"]
    ).reset_index()
    fig = plt.figure(figsize=(30, 15))
    for i, study in enumerate(rna_bioprojects):
        ax = fig.add_subplot(1, 3, i + 1)
        sns.boxenplot(
            ax=ax,
            data=data[data.study == study],
            x="ra_at_1in1000",
            y="pathogen",
            order=plotting_order.pathogen.unique(),
            hue="location",
            showfliers=False,
        )
        for num_reads, artist_list in zip(
            plotting_order[plotting_order.study == study].viral_reads,
            ax.collections,
        ):
            artist_list.set_alpha(min(num_reads / 10 + 0.02, 1.0))
        ax.set_xscale("log")
        ax.set_ylabel("")
    return fig


def count_viral_reads(
    df: pd.DataFrame, by_location: bool = False
) -> pd.DataFrame:
    groups = ["pathogen", "predictor_type", "study"]
    if by_location:
        groups.append("fine_location")
    out = df.groupby(groups).viral_reads.sum().reset_index()
    out["reads_by_pathogen"] = out.viral_reads.groupby(out.pathogen).transform(
        "sum"
    )
    return out


def save_plot(fig, figdir: Path, name: str) -> None:
    for ext in ["pdf", "png"]:
        fig.savefig(figdir / f"{name}.{ext}", bbox_inches="tight")


def start() -> None:
    figdir = Path("fig")
    figdir.mkdir(exist_ok=True)
    fits_df = pd.read_csv("fits.tsv", sep="\t")
    input_df = pd.read_csv("input.tsv", sep="\t")
    fig_overall = plot_overall(fits_df, input_df)
    save_plot(fig_overall, figdir, "overall-boxen")
    fig_by_location = plot_by_location(fits_df, input_df)
    save_plot(fig_by_location, figdir, "by_location-boxen")


if __name__ == "__main__":
    start()
