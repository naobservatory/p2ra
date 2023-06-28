#!/usr/bin/env python3
from pathlib import Path

import matplotlib.pyplot as plt  # type: ignore
import pandas as pd
import seaborn as sns  # type: ignore

from mgs import rna_bioprojects

study_order = ["Crits-Christoph", "Rothman", "Spurbeck"]

FIGSIZE = (4, 6)


def separate_viruses(ax) -> None:
    yticks = ax.get_yticks()
    ax.hlines(
        [(y1 + y2) / 2 for y1, y2 in zip(yticks[:-1], yticks[1:])],
        *ax.get_xlim(),
        linestyle="dashed",
        color="0.5",
        linewidth=0.5,
    )


def adjust_axes(ax) -> None:
    yticks = ax.get_yticks()
    # Y-axis is reflected
    ax.set_ylim([max(yticks) + 0.5, min(yticks - 0.5)])
    ax.tick_params(left=False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.set_xscale("log")
    ax.set_xlabel(r"RA$_{1:1000}$")
    ax.set_ylabel("")
    ax.legend(
        title="MGS study",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        borderaxespad=0,
        frameon=False,
    )


def plot_overall(data: pd.DataFrame, input_data: pd.DataFrame) -> plt.Figure:
    viral_reads = count_viral_reads(input_data)
    plotting_order = viral_reads.sort_values(
        ["reads_by_pathogen", "pathogen"]
    ).reset_index()
    fig = plt.figure(figsize=FIGSIZE)
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
    separate_viruses(ax)
    adjust_axes(ax)
    return fig


def plot_virus(
    data: pd.DataFrame, input_data: pd.DataFrame, pathogen: str
) -> plt.Figure:
    viral_reads = count_viral_reads(
        input_data[input_data.pathogen == pathogen], by_location=True
    )
    fig = plt.figure(figsize=FIGSIZE)
    ax = fig.add_subplot()
    sns.boxenplot(
        ax=ax,
        data=data[(data.location != "Overall") & (data.pathogen == pathogen)],
        x="ra_at_1in1000",
        y="location",
        hue="study",
        showfliers=False,
    )
    for num_reads, artist_list in zip(
        viral_reads.viral_reads,
        ax.collections,
    ):
        artist_list.set_alpha(min(num_reads / 10 + 0.02, 1.0))
    adjust_axes(ax)
    ax.set_title(pathogen)
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
    fits_df["study"] = [x.replace("_", "-").title() for x in fits_df["study"]]
    input_df = pd.read_csv("input.tsv", sep="\t")
    fig_overall = plot_overall(fits_df, input_df)
    save_plot(fig_overall, figdir, "overall-boxen")
    for pathogen in input_df.pathogen.unique():
        fig_virus = plot_virus(fits_df, input_df, pathogen)
        save_plot(fig_virus, figdir, f"{pathogen.replace(' ', '_')}-boxen")


if __name__ == "__main__":
    start()
