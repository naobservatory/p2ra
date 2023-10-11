#!/usr/bin/env python3
from pathlib import Path
import sys

sys.path.append("..")


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
    ax.vlines(
        ax.get_xticks()[1:-1],
        *ax.get_ylim(),
        linewidth=0.05,
        color="black",
        zorder=-1,
    )
    # ax.set_xscale("log")
    ax.set_xlabel(
        r"$\mathrm{RA}"
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


def plot_three_virus(
    data: pd.DataFrame,
    input_data: pd.DataFrame,
    viruses: dict[str, tuple[float, float]],
    predictor_type: str,
    fig_title: str,
    axes: list[plt.Axes, plt.Axes, plt.Axes],
) -> list[plt.Axes, plt.Axes, plt.Axes]:
    final_axes = []
    for i, ((pathogen, xlim), ax) in enumerate(zip(viruses.items(), axes)):
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
        num_spurbeck = 10
        num_rothman = 8
        num_crits_christoph = 4
        ax.hlines(
            [
                num_spurbeck - 0.5,
                num_spurbeck + num_rothman - 0.5,
                num_spurbeck + num_rothman + num_crits_christoph - 0.5,
            ],
            *ax.get_xlim(),
            linestyle="solid",
            color="k",
            linewidth=0.5,
        )
        if i == 0:
            sub_fig_title = ax.set_title(
                fig_title,
                fontweight="bold",
                loc="left",
            )
            sub_fig_title.set_position((-0.35, 1.05))
        if i == 2:
            x_text = ax.get_xlim()[1] + 0.1
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
            if predictor_type == "prevalence":
                ax.text(
                    x_text,
                    num_spurbeck + num_rothman + num_crits_christoph - 0.4,
                    "Brinch",
                    va="top",
                )

        adjust_axes(ax, predictor_type=predictor_type)
        plot_title = ax.set_title(pathogen)
        ax.get_legend().remove()
        if i != 1:
            ax.set_xlabel("")
        if i > 0:
            ax.set_yticklabels([])

        final_axes.append(ax)
    return final_axes


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


def composite_figure(
    data: pd.DataFrame,
    input_data: pd.DataFrame,
) -> plt.Figure:
    fig = plt.figure(
        figsize=(8, 12),
    )

    gs = fig.add_gridspec(2, 3, height_ratios=[7, 7], hspace=0.2)

    incidence_viruses = {
        "Norovirus (GII)": (-9.0, -2.0),
        "Norovirus (GI)": (-9.0, -2.0),
        "SARS-COV-2": (-11.0, -5.0),
    }
    plot_three_virus(
        data,
        input_data,
        incidence_viruses,
        "incidence",
        "a",
        [
            fig.add_subplot(gs[0, 0]),
            fig.add_subplot(gs[0, 1]),
            fig.add_subplot(gs[0, 2]),
        ],
    )

    prevalence_viruses = {
        "JCV": (-14.0, -7.0),
        "BKV": (-14.0, -7.0),
        "MCV": (-14.0, -7.0),
    }
    plot_three_virus(
        data,
        input_data,
        prevalence_viruses,
        "prevalence",
        "b",
        [
            fig.add_subplot(gs[1, 0]),
            fig.add_subplot(gs[1, 1]),
            fig.add_subplot(gs[1, 2]),
        ],
    )
    return fig


def save_plot(fig, figdir: Path, name: str) -> None:
    for ext in ["pdf", "png"]:
        fig.savefig(figdir / f"{name}.{ext}", bbox_inches="tight", dpi=600)


def start() -> None:
    parent_dir = Path("..")
    figdir = Path(parent_dir / "figures")
    figdir.mkdir(exist_ok=True)
    fits_df = pd.read_csv(parent_dir / "fits.tsv", sep="\t")
    fits_df["study"] = fits_df.study.map(study_name)
    fits_df["log10ra"] = np.log10(fits_df.ra_at_1in100)
    input_df = pd.read_csv(parent_dir / "input.tsv", sep="\t")
    input_df["study"] = input_df.study.map(study_name)
    # TODO: Store these in the files instead?
    input_df["nucleic_acid"] = input_df.pathogen.map(nucleic_acid)
    input_df["selection_round"] = input_df.pathogen.map(selection_round)
    input_df["observed?"] = input_df.viral_reads > 0
    # For consistency between dataframes (TODO: fix that elsewhere)
    input_df["location"] = input_df.fine_location

    fig = composite_figure(fits_df, input_df)

    save_plot(fig, figdir, "composite_fig_5")


if __name__ == "__main__":
    start()
