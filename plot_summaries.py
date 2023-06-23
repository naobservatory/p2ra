#!/usr/bin/env python3
from pathlib import Path

import matplotlib.pyplot as plt  # type: ignore
import pandas as pd
import seaborn as sns  # type: ignore

from mgs import rna_bioprojects


def plot_overall(data: pd.DataFrame) -> plt.Figure:
    fig = plt.figure(figsize=(6, 10))
    ax = fig.add_subplot()
    sns.boxenplot(
        ax=ax,
        data=data[data.location == "Overall"],
        x="ra_at_1in1000",
        y="path_tax",
        hue="study",
        showfliers=False,
    )
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


def start() -> None:
    figdir = Path("fig")
    figdir.mkdir(exist_ok=True)
    datafile = "fits.tsv"
    df = pd.read_csv(datafile, sep="\t")
    df["path_tax"] = df.pathogen + "-" + df.taxids.astype("string")
    fig_overall = plot_overall(df)
    fig_overall.savefig(figdir / "overall-boxen.pdf", bbox_inches="tight")
    fig_overall.savefig(figdir / "overall-boxen.png", bbox_inches="tight")
    fig_by_location = plot_by_location(df)
    fig_by_location.savefig(
        figdir / "by_location-boxen.pdf", bbox_inches="tight"
    )
    fig_by_location.savefig(
        figdir / "by_location-boxen.png", bbox_inches="tight"
    )


if __name__ == "__main__":
    start()
