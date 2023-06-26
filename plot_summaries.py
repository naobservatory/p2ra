#!/usr/bin/env python3
from pathlib import Path

import matplotlib.pyplot as plt  # type: ignore
import pandas as pd
import seaborn as sns  # type: ignore

from mgs import rna_bioprojects


def plot_overall(data: pd.DataFrame, input_data: pd.DataFrame) -> plt.Figure:
    viral_reads = count_viral_reads(input_data)
    print(len(viral_reads))
    print(sorted(viral_reads))
    fig = plt.figure(figsize=(6, 10))
    ax = fig.add_subplot()
    sns.boxenplot(
        ax=ax,
        data=data[data.location == "Overall"],
        x="ra_at_1in1000",
        y="path_tax",
        order=sorted(data.path_tax.unique()),
        hue="study",
        showfliers=False,
    )
    for num_reads, artist_list in zip(viral_reads, ax.collections):
        artist_list.set_alpha(min(num_reads / 10, 1.0))
    ax.set_xscale("log")
    return fig


def plot_by_location(data: pd.DataFrame) -> plt.Figure:
    fig = plt.figure(figsize=(30, 15))
    for i, study in enumerate(rna_bioprojects):
        ax = fig.add_subplot(1, 3, i + 1)
        sns.boxenplot(
            ax=ax,
            data=data[data.study == study],
            x="ra_at_1in1000",
            y="path_tax",
            hue="location",
            showfliers=False,
        )
        ax.set_xscale("log")
    return fig


def count_viral_reads(
    df: pd.DataFrame, by_location: bool = False
) -> pd.Series:
    groups = ["pathogen", "taxids", "predictor_type", "study"]
    if by_location:
        groups.append("fine_location")
    return df.groupby(groups).viral_reads.sum()


def save_plot(fig, figdir: Path, name: str) -> None:
    for ext in ["pdf", "png"]:
        fig.savefig(figdir / f"{name}.{ext}", bbox_inches="tight")


def start() -> None:
    figdir = Path("fig")
    figdir.mkdir(exist_ok=True)
    fits_df = pd.read_csv("fits.tsv", sep="\t")
    input_df = pd.read_csv("input.tsv", sep="\t")
    print(count_viral_reads(input_df))
    fits_df["path_tax"] = (
        fits_df.pathogen + "-" + fits_df.taxids.astype("string")
    )
    input_df["path_tax"] = (
        input_df.pathogen + "-" + input_df.taxids.astype("string")
    )
    fig_overall = plot_overall(fits_df, input_df)
    save_plot(fig_overall, figdir, "overall-boxen")
    fig_by_location = plot_by_location(fits_df)
    save_plot(fig_by_location, figdir, "by_location-boxen")


if __name__ == "__main__":
    start()
