#!/usr/bin/env python3
from pathlib import Path

import matplotlib.pyplot as plt  # type: ignore
import pandas as pd
import seaborn as sns  # type: ignore

from pathogens import pathogens


def nucleic_acid(pathogen: str) -> str:
    return pathogens[pathogen].pathogen_chars.na_type.value


def selection_round(pathogen: str) -> str:
    return pathogens[pathogen].pathogen_chars.selection.value


def study_name(study: str) -> str:
    return {
        "brinch": "Brinch (DNA)",
        "crits_christoph": "Crits-Christoph",
        "rothman": "Rothman",
        "spurbeck": "Spurbeck",
    }[study]


plt.rcParams["font.size"] = 8


def separate_viruses(ax) -> None:
    yticks = ax.get_yticks()
    ax.hlines(
        [(y1 + y2) / 2 for y1, y2 in zip(yticks[:-1], yticks[1:])],
        *ax.get_xlim(),
        linestyle="dotted",
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
        "Expected relative abundance when "
        f"{predictor_type} = 1:1000, "
        r"$RA"
        f"{predictor_type[0]}"
        r"(1\perthousand)$"
    )
    ax.set_ylabel("")


def plot_boxen(
    ax,
    data: pd.DataFrame,
    viral_reads: pd.DataFrame,
    y: str,
    sorting_order: list[str],
    ascending: list[bool],
    show_zero_counts: bool = False,
) -> None:
    assert len(sorting_order) == len(ascending)
    plotting_order = viral_reads.sort_values(
        sorting_order, ascending=ascending
    ).reset_index()
    sns.boxenplot(
        ax=ax,
        data=data,
        x="ra_at_1in1000",
        y=y,
        order=plotting_order[y].unique(),
        hue="study",
        hue_order=plotting_order.study.unique(),
        showfliers=False,
        box_kws={"linewidth": 0.5},
    )

    lhs = ax.get_xlim()[0]
    for num_reads, line, patches in zip(
        plotting_order.viral_reads, ax.lines, ax.collections
    ):
        alpha = min((num_reads + 1) / 10, 1.0)
        line.set_alpha(alpha)
        patches.set_alpha(alpha)
        if (not show_zero_counts) and (num_reads == 0):
            patches.set_visible(False)
            line.set_visible(False)
            # The third patch from center marks the ~94 percentile
            path = patches.get_paths()[-3]
            right_edge_midpoint = (
                (path.vertices[1][0] + path.vertices[2][0]) / 2,
                (path.vertices[1][1] + path.vertices[2][1]) / 2,
            )
            color = patches.get_cmap()(0)
            ax.hlines(
                right_edge_midpoint[1],
                lhs,
                right_edge_midpoint[0],
                linestyle="dashed",
                color=color,
                alpha=0.1,
            )
            ax.scatter(
                *right_edge_midpoint,
                color=color,
                alpha=0.5,
                marker="<",
                s=10,
            )


def plot_incidence(data: pd.DataFrame, input_data: pd.DataFrame) -> plt.Figure:
    predictor_type = "incidence"
    fig = plt.figure(figsize=(4, 3))
    ax = fig.add_subplot(1, 1, 1)
    plot_boxen(
        ax=ax,
        data=data[
            (data.predictor_type == predictor_type)
            & (data.location == "Overall")
        ],
        viral_reads=count_viral_reads(
            input_data[input_data.predictor_type == predictor_type]
        ),
        y="tidy_name",
        sorting_order=[
            "nucleic_acid",
            "selection_round",
            "samples_observed_by_tidy_name",
            "tidy_name",
            "study",
        ],
        ascending=[False, True, False, True, False],
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
    fig = plt.figure(figsize=(4, 6))
    ax = fig.add_subplot(1, 1, 1)
    plot_boxen(
        ax=ax,
        data=data[
            (data.predictor_type == predictor_type)
            & (data.location == "Overall")
        ],
        viral_reads=count_viral_reads(
            input_data[input_data.predictor_type == predictor_type]
        ),
        y="tidy_name",
        sorting_order=[
            "nucleic_acid",
            "selection_round",
            "samples_observed_by_tidy_name",
            "tidy_name",
            "study",
        ],
        ascending=[False, True, False, True, False],
    )
    ax.set_xlim([1e-16, 1e-8])
    separate_viruses(ax)
    # TODO Get these values automatically
    num_rna_1 = 2
    num_dna_1 = 4
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
    data: pd.DataFrame,
    input_data: pd.DataFrame,
    viruses: dict[str, tuple[float, float]],
    predictor_type: str,
) -> plt.Figure:
    fig = plt.figure(figsize=(6, 4))
    for i, (pathogen, xlim) in enumerate(viruses.items()):
        ax = fig.add_subplot(1, 3, i + 1)
        plot_boxen(
            ax=ax,
            data=data[
                (data.location != "Overall") & (data.tidy_name == pathogen)
            ],
            viral_reads=count_viral_reads(
                input_data[input_data.tidy_name == pathogen], by_location=True
            ),
            y="location",
            sorting_order=["study", "location"],
            ascending=[False, True],
            show_zero_counts=True,
        )
        ax.set_xlim(xlim)
        # TODO Get these values automatically
        num_spurbeck = 10
        num_rothman = 8
        ax.hlines(
            [num_spurbeck - 0.5, num_spurbeck + num_rothman - 0.5],
            *ax.get_xlim(),
            linestyle="solid",
            color="k",
            linewidth=0.5,
        )
        if i == 2:
            x_text = 1.2e-6
            ax.text(x_text, -0.4, "Spurbeck", va="top")
            ax.text(
                x_text,
                num_spurbeck - 0.4,
                "Rothman",
                va="top",
            )
            ax.text(
                x_text,
                num_spurbeck + num_rothman - 0.4,
                "Crits-Christoph",
                va="top",
            )
        adjust_axes(ax, predictor_type=predictor_type)
        ax.set_title(pathogen)
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
        groups.append("location")
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
    fits_df["study"] = fits_df.study.map(study_name)
    input_df = pd.read_csv("input.tsv", sep="\t")
    input_df["study"] = input_df.study.map(study_name)
    # TODO: Store these in the files instead?
    input_df["nucleic_acid"] = input_df.pathogen.map(nucleic_acid)
    input_df["selection_round"] = input_df.pathogen.map(selection_round)
    input_df["observed?"] = input_df.viral_reads > 0
    # For consistency between dataframes (TODO: fix that elsewhere)
    input_df["location"] = input_df.fine_location
    fig_incidence = plot_incidence(fits_df, input_df)
    save_plot(fig_incidence, figdir, "incidence-boxen")
    fig_prevalence = plot_prevalence(fits_df, input_df)
    save_plot(fig_prevalence, figdir, "prevalence-boxen")
    incidence_viruses = {
        "Norovirus (GII)": (1e-10, 1e-3),
        "Norovirus (GI)": (1e-10, 1e-3),
        "SARS-COV-2": (1e-12, 1e-6),
    }
    fig_three_virus_incidence = plot_three_virus(
        fits_df, input_df, incidence_viruses, "incidence"
    )
    save_plot(fig_three_virus_incidence, figdir, "by_location_incidence-boxen")


if __name__ == "__main__":
    start()
