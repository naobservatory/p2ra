#!/usr/bin/env python3
from pathlib import Path

from pathogens import pathogens

import matplotlib.pyplot as plt  # type: ignore
import pandas as pd
import seaborn as sns  # type: ignore

nucleic_acid = {
    pathogen_name: module.pathogen_chars.na_type.value
    for pathogen_name, module in pathogens.items()
}


def selection_round(pathogen: str) -> int:
    if pathogen in ["jcv", "bkv", "mcv", "aav2", "aav5", "aav6"]:
        return 2
    else:
        return 1


plt.rcParams["font.size"] = 8

FIGSIZE = (4, 6)


def separate_viruses(ax) -> None:
    yticks = ax.get_yticks()
    ax.hlines(
        [(y1 + y2) / 2 for y1, y2 in zip(yticks[:-1], yticks[1:])],
        *ax.get_xlim(),
        linestyle="dotted",
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


def adjust_axes(ax, predictor_type: str) -> None:
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
        f"{predictor_type} = 1:1000 "
        r"($RA"
        f"^{predictor_type[0]}"
        r"_{1:1000}$)"
    )
    ax.set_ylabel("")


def plot_boxen(ax, data: pd.DataFrame, input_data: pd.DataFrame) -> None:
    viral_reads = count_viral_reads(input_data)
    print(viral_reads)
    plotting_order = viral_reads.sort_values(
        [
            "nucleic_acid",
            "selection_round",
            "samples_observed_by_tidy_name",
            "tidy_name",
            "study",
        ],
        ascending=[False, True, False, True, False],
    ).reset_index()
    sns.boxenplot(
        ax=ax,
        data=data[data.location == "Overall"],
        x="ra_at_1in1000",
        y="tidy_name",
        order=plotting_order.tidy_name.unique(),
        hue="study",
        hue_order=plotting_order.study.unique(),
        showfliers=False,
        box_kws={"linewidth": 0.5},
    )
    for num_reads, artist_list in zip(
        plotting_order.viral_reads, ax.collections
    ):
        artist_list.set_alpha(min(num_reads / 10 + 0.02, 1.0))


def plot_incidence(data: pd.DataFrame, input_data: pd.DataFrame) -> plt.Figure:
    predictor_type = "incidence"
    fig = plt.figure(figsize=(4, 3))
    ax = fig.add_subplot(1, 1, 1)
    plot_boxen(
        ax,
        data[data.predictor_type == predictor_type],
        input_data[input_data.predictor_type == predictor_type],
    )
    ax.set_xlim([1e-14, 1e-4])
    separate_viruses(ax)
    adjust_axes(ax, predictor_type=predictor_type)
    ax.legend(
        title="MGS study",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        borderaxespad=0,
        frameon=False,
    )
    return fig


def plot_prevalence(
    data: pd.DataFrame, input_data: pd.DataFrame
) -> plt.Figure:
    predictor_type = "prevalence"
    fig = plt.figure(figsize=FIGSIZE)
    ax = fig.add_subplot(1, 1, 1)
    plot_boxen(
        ax,
        data[data.predictor_type == predictor_type],
        input_data[input_data.predictor_type == predictor_type],
    )
    ax.set_xlim([1e-16, 1e-8])
    separate_viruses(ax)
    # TODO Get these values automatically
    num_rna_1 = 2
    num_dna_1 = 6
    ax.hlines(
        [num_rna_1 - 0.5, num_rna_1 + num_dna_1 - 0.5],
        *ax.get_xlim(),
        linestyle="solid",
        color="k",
        linewidth=0.5,
    )
    ax.text(1.1e-8, -0.4, "RNA viruses\nSelection Round 1", va="top")
    ax.text(
        1.1e-8, num_rna_1 - 0.4, "DNA viruses\nSelection Round 1", va="top"
    )
    ax.text(
        1.1e-8,
        num_rna_1 + num_dna_1 - 0.4,
        "DNA viruses\nSelection Round 2",
        va="top",
    )
    adjust_axes(ax, predictor_type=predictor_type)
    ax.legend(
        title="MGS study",
        bbox_to_anchor=(1.02, 0),
        loc="lower left",
        borderaxespad=0,
        frameon=False,
    )
    return fig


def plot_three_virus(
    data: pd.DataFrame, input_data: pd.DataFrame
) -> plt.Figure:
    viruses = {
        "Norovirus (GII)": (1e-10, 1e-3),
        "Norovirus (GI)": (1e-10, 1e-3),
        "SARS-COV-2": (1e-13, 1e-5),
    }
    fig = plt.figure(figsize=(6, 4))
    for i, (pathogen, xlim) in enumerate(viruses.items()):
        viral_reads = count_viral_reads(
            input_data[input_data.tidy_name == pathogen], by_location=True
        )
        ax = fig.add_subplot(1, 3, i + 1)
        to_plot = data[
            (data.location != "Overall") & (data.tidy_name == pathogen)
        ]
        sns.boxenplot(
            ax=ax,
            data=to_plot,
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
        adjust_axes(ax, predictor_type="incidence")
        ax.set_title(pathogen)
        if i == 2:
            ax.legend(
                title="MGS study",
                bbox_to_anchor=(1.02, 1),
                loc="upper left",
                borderaxespad=0,
                frameon=False,
            )
        else:
            ax.get_legend().remove()
        if i != 1:
            ax.set_xlabel("")
        if i > 0:
            ax.set_yticklabels([])
    return fig


def count_viral_reads(
    df: pd.DataFrame, by_location: bool = False
) -> pd.DataFrame:
    groups = [
        "pathogen",
        "tidy_name",
        "predictor_type",
        "study",
        "nucleic_acid",
        "selection_round",
    ]
    if by_location:
        groups.append("fine_location")
    out = df.groupby(groups)[["viral_reads", "observed?"]].sum().reset_index()
    out["reads_by_tidy_name"] = out.viral_reads.groupby(
        out.tidy_name
    ).transform("sum")
    out["samples_observed_by_tidy_name"] = (
        out["observed?"].groupby(out.tidy_name).transform("sum")
    )
    return out


def save_plot(fig, figdir: Path, name: str) -> None:
    for ext in ["pdf", "png"]:
        fig.savefig(figdir / f"{name}.{ext}", bbox_inches="tight")


def start() -> None:
    figdir = Path("fig")
    figdir.mkdir(exist_ok=True)
    fits_df = pd.read_csv("fits.tsv", sep="\t")
    studies = {
        "brinch": "Brinch (DNA)",
        "crits_christoph": "Crits-Christoph",
        "rothman": "Rothman",
        "spurbeck": "Spurbeck",
    }
    fits_df["study"] = fits_df.study.map(lambda s: studies[s])
    input_df = pd.read_csv("input.tsv", sep="\t")
    input_df["study"] = input_df.study.map(lambda s: studies[s])
    input_df["nucleic_acid"] = input_df.pathogen.map(lambda p: nucleic_acid[p])
    input_df["selection_round"] = input_df.pathogen.map(selection_round)
    input_df["observed?"] = input_df.viral_reads > 0
    print(input_df)
    fig_incidence = plot_incidence(fits_df, input_df)
    save_plot(fig_incidence, figdir, "incidence-boxen")
    fig_prevalence = plot_prevalence(fits_df, input_df)
    save_plot(fig_prevalence, figdir, "prevalence-boxen")
    fig_three_virus = plot_three_virus(fits_df, input_df)
    save_plot(fig_three_virus, figdir, "three_virus-boxen")


if __name__ == "__main__":
    start()
