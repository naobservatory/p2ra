#!/usr/bin/env python3
from pathlib import Path

import matplotlib.pyplot as plt  # type: ignore
import pandas as pd
import seaborn as sns  # type: ignore

study_order = ["Crits-Christoph", "Rothman", "Spurbeck", "Brinch"]

plt.rcParams["font.size"] = 8

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


def separate_studies(ax) -> None:
    # yticks = ax.get_yticks()
    ax.hlines(
        [3.5, 11.5],  # TODO: get this from data
        # [(y1 + y2) / 2 for y1, y2 in zip(yticks[:-1], yticks[1:])],
        *ax.get_xlim(),
        linestyle="dashed",
        color="0.5",
        linewidth=0.5,
    )


def italicize_incidence_viruses(ax) -> None:
    # TODO: Get from data
    incidence_viruses = ["SARS-COV-2", "Norovirus (GI)", "Norovirus (GII)"]
    for label in ax.get_yticklabels():
        if label.get_text() in incidence_viruses:
            label.set_fontstyle("oblique")


def adjust_axes(ax) -> None:
    yticks = ax.get_yticks()
    # Y-axis is reflected
    ax.set_ylim([max(yticks) + 0.5, min(yticks - 0.5)])
    ax.tick_params(left=False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.set_xscale("log")
    ax.set_xlabel(
        "Expected relative abundance when\n"
        "incidence|prevalence = 1:1000 "
        r"($RA_{1:1000}$)"
    )
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
        ["reads_by_pathogen", "tidy_name"]
    ).reset_index()
    fig = plt.figure(figsize=FIGSIZE)
    ax = fig.add_subplot(1, 1, 1)
    sns.boxenplot(
        ax=ax,
        data=data[data.location == "Overall"],
        x="ra_at_1in1000",
        y="tidy_name",
        order=plotting_order.tidy_name.unique(),
        hue="study",
        hue_order=study_order,
        showfliers=False,
        box_kws={"linewidth": 0.5},
    )
    for num_reads, artist_list in zip(
        plotting_order.viral_reads, ax.collections
    ):
        artist_list.set_alpha(min(num_reads / 10 + 0.02, 1.0))
    ax.set_xlim([1e-18, 1e-3])
    separate_viruses(ax)
    adjust_axes(ax)
    italicize_incidence_viruses(ax)
    return fig


def plot_virus(
    data: pd.DataFrame, input_data: pd.DataFrame, pathogen: str
) -> plt.Figure:
    viral_reads = count_viral_reads(
        input_data[input_data.pathogen == pathogen], by_location=True
    )
    fig = plt.figure(figsize=(4, 3))
    ax = fig.add_subplot()
    sns.boxenplot(
        ax=ax,
        data=data[(data.location != "Overall") & (data.pathogen == pathogen)],
        x="ra_at_1in1000",
        y="location",
        hue="study",
        showfliers=False,
        box_kws={"linewidth": 0.5},
    )
    for num_reads, artist_list in zip(
        viral_reads.viral_reads,
        ax.collections,
    ):
        artist_list.set_alpha(min(num_reads / 10 + 0.02, 1.0))
    separate_studies(ax)
    adjust_axes(ax)
    ax.set_title(pathogen)
    return fig


def plot_three_virus(
    data: pd.DataFrame, input_data: pd.DataFrame
) -> plt.Figure:
    viruses = {
        "Norovirus (GII)": (1e-10, 1e-3),
        "MCV": (1e-16, 1e-9),
        "SARS-COV-2": (1e-13, 1e-5),
    }
    fig = plt.figure(figsize=(6, 4))
    for i, (pathogen, xlim) in enumerate(viruses.items()):
        viral_reads = count_viral_reads(
            input_data[input_data.pathogen == pathogen], by_location=True
        )
        ax = fig.add_subplot(1, 3, i + 1)
        sns.boxenplot(
            ax=ax,
            data=data[
                (data.location != "Overall") & (data.pathogen == pathogen)
            ],
            x="ra_at_1in1000",
            y="location",
            hue="study",
            showfliers=False,
            box_kws={"linewidth": 0.5},
        )
        for num_reads, artist_list in zip(
            viral_reads.viral_reads,
            ax.collections,
        ):
            artist_list.set_alpha(min(num_reads / 10 + 0.02, 1.0))
        ax.set_xlim(xlim)
        separate_studies(ax)
        adjust_axes(ax)
        ax.set_title(pathogen)
        if i < 2:
            ax.get_legend().remove()
        if i != 1:
            ax.set_xlabel("")
        if i > 0:
            ax.set_yticklabels([])
    return fig


def count_viral_reads(
    df: pd.DataFrame, by_location: bool = False
) -> pd.DataFrame:
    groups = ["pathogen", "tidy_name", "predictor_type", "study"]
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
    fig_three_virus = plot_three_virus(fits_df, input_df)
    save_plot(fig_three_virus, figdir, "three_virus-boxen")
    # for pathogen in input_df.pathogen.unique():
    #     fig_virus = plot_virus(fits_df, input_df, pathogen)
    #     save_plot(fig_virus, figdir, f"{pathogen.replace(' ', '_')}-boxen")


if __name__ == "__main__":
    start()
