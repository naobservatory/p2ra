#!/usr/bin/env python3
from pathlib import Path

import matplotlib.patches as mpatches  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
import matplotlib.ticker as ticker  # type: ignore
import numpy as np
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
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(format_func))
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    # ax.set_xscale("log")
    ax.set_xlabel(
        r"$RA"
        f"{predictor_type[0]}"
        r"(1\%)$"
        ": expected relative abundance at 1% "
        f"{predictor_type} "
    )
    ax.set_ylabel("")


def plot_violin(
    ax,
    data: pd.DataFrame,
    viral_reads: pd.DataFrame,
    y: str,
    sorting_order: list[str],
    ascending: list[bool],
    hatch_zero_counts: bool = False,
    violin_scale=1.0,
) -> None:
    assert len(sorting_order) == len(ascending)
    plotting_order = viral_reads.sort_values(
        sorting_order, ascending=ascending
    ).reset_index()
    sns.violinplot(
        ax=ax,
        data=data,
        x="log10ra",
        y=y,
        order=plotting_order[y].unique(),
        hue="study",
        hue_order=plotting_order.study.unique(),
        inner=None,
        linewidth=0.0,
        bw=0.5,
        scale="area",
        scale_hue=False,
        cut=0,
    )
    x_min = ax.get_xlim()[0]
    for num_reads, patches in zip(plotting_order.viral_reads, ax.collections):
        # alpha = min((num_reads + 1) / 10, 1.0)
        if num_reads == 0:
            alpha = 0.5
        elif num_reads < 10:
            alpha = 0.5
        else:
            alpha = 1.0
        patches.set_alpha(alpha)
        # Make violins fatter and hatch if zero counts
        for path in patches.get_paths():
            y_mid = path.vertices[0, 1]
            path.vertices[:, 1] = (
                violin_scale * (path.vertices[:, 1] - y_mid) + y_mid
            )
            if (not hatch_zero_counts) and (num_reads == 0):
                color = patches.get_facecolor()
                y_max = np.max(path.vertices[:, 1])
                y_min = np.min(path.vertices[:, 1])
                x_max = path.vertices[np.argmax(path.vertices[:, 1]), 0]
                rect = mpatches.Rectangle(
                    (x_min, y_min),
                    x_max - x_min,
                    y_max - y_min,
                    facecolor=color,
                    linewidth=0.0,
                    alpha=alpha,
                    fill=False,
                    hatch="|||",
                    edgecolor=color,
                )
                ax.add_patch(rect)


def format_func(value, tick_number):
    return r"$10^{{{}}}$".format(int(value))


def plot_incidence(data: pd.DataFrame, input_data: pd.DataFrame) -> plt.Figure:
    predictor_type = "incidence"
    fig = plt.figure(figsize=(4, 3))
    ax = fig.add_subplot(1, 1, 1)

    plot_violin(
        ax=ax,
        data=data[
            (data.predictor_type == predictor_type)
            & (data.location == "Overall")
            & ~(
                (data.study == "Crits-Christoph")
                & (data.pathogen == "influenza")
            )
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
        violin_scale=2.0,
    )
    ax.set_xlim([-13, -3])
    separate_viruses(ax)
    adjust_axes(ax, predictor_type=predictor_type)
    legend = ax.legend(
        title="MGS study",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        borderaxespad=0,
        frameon=False,
    )
    for legend_handle in legend.legend_handles:
        legend_handle.set_edgecolor(legend_handle.get_facecolor())
    return fig


def plot_prevalence(
    data: pd.DataFrame, input_data: pd.DataFrame
) -> plt.Figure:
    predictor_type = "prevalence"
    fig = plt.figure(figsize=(4, 6))
    ax = fig.add_subplot(1, 1, 1)
    plot_violin(
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
        violin_scale=1.5,
    )
    ax.set_xlim([-15, -7])
    ax.set_xticks(list(range(-15, -5, 2)))
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
    text_x = np.log10(1.1e-7)
    ax.text(text_x, -0.4, "RNA viruses\nSelection Round 1", va="top")
    ax.text(
        text_x, num_rna_1 - 0.4, "DNA viruses\nSelection Round 1", va="top"
    )
    ax.text(
        text_x,
        num_rna_1 + num_dna_1 - 0.4,
        "DNA viruses\nSelection Round 2",
        va="top",
    )
    adjust_axes(ax, predictor_type=predictor_type)
    legend = ax.legend(
        title="MGS study",
        bbox_to_anchor=(1.02, 0),
        loc="lower left",
        borderaxespad=0,
        frameon=False,
    )
    for legend_handle in legend.legend_handles:
        legend_handle.set_edgecolor(legend_handle.get_facecolor())
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
        plot_violin(
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
            violin_scale=2.5,
            hatch_zero_counts=True,
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
            x_text = np.log10(1.2e-5)
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
        fig.savefig(figdir / f"{name}.{ext}", bbox_inches="tight", dpi=600)


def start() -> None:
    figdir = Path("fig")
    figdir.mkdir(exist_ok=True)
    fits_df = pd.read_csv("fits.tsv", sep="\t")
    fits_df["study"] = fits_df.study.map(study_name)
    fits_df["log10ra"] = np.log10(fits_df.ra_at_1in100)
    input_df = pd.read_csv("input.tsv", sep="\t")
    input_df["study"] = input_df.study.map(study_name)
    # TODO: Store these in the files instead?
    input_df["nucleic_acid"] = input_df.pathogen.map(nucleic_acid)
    input_df["selection_round"] = input_df.pathogen.map(selection_round)
    input_df["observed?"] = input_df.viral_reads > 0
    # For consistency between dataframes (TODO: fix that elsewhere)
    input_df["location"] = input_df.fine_location
    fig_incidence = plot_incidence(fits_df, input_df)
    save_plot(fig_incidence, figdir, "incidence-violin")
    fig_prevalence = plot_prevalence(fits_df, input_df)
    save_plot(fig_prevalence, figdir, "prevalence-violin")
    incidence_viruses = {
        "Norovirus (GII)": (-9.0, -2.0),
        "Norovirus (GI)": (-9.0, -2.0),
        "SARS-COV-2": (-11.0, -5.0),
    }
    fig_three_virus_incidence = plot_three_virus(
        fits_df, input_df, incidence_viruses, "incidence"
    )
    save_plot(
        fig_three_virus_incidence, figdir, "by_location_incidence-violin"
    )


if __name__ == "__main__":
    start()
