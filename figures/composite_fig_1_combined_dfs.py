import os
import gzip
import json
import subprocess
import matplotlib.pyplot as plt  # type: ignore
from matplotlib.gridspec import GridSpec
from matplotlib.colors import ListedColormap
import matplotlib.cbook as cbook
import pandas as pd
from scipy.stats import gmean
import numpy as np
from PIL import Image
import seaborn as sns


dashboard = os.path.expanduser("~/code/mgs-pipeline/dashboard/")


with open(os.path.join(dashboard, "human_virus_sample_counts.json")) as inf:
    human_virus_sample_counts = json.load(inf)

with open(os.path.join(dashboard, "metadata_samples.json")) as inf:
    metadata_samples = json.load(inf)

with open(os.path.join(dashboard, "metadata_bioprojects.json")) as inf:
    metadata_bioprojects = json.load(inf)

with open(os.path.join(dashboard, "metadata_papers.json")) as inf:
    metadata_papers = json.load(inf)

with open(os.path.join(dashboard, "taxonomic_names.json")) as inf:
    taxonomic_names = json.load(inf)


studies = list(metadata_papers.keys())


target_taxa = {
    10239: ("viruses", "Viruses"),
    2731341: ("duplodnaviria", "DNA Viruses"),
    2732004: ("varidnaviria", "DNA Viruses"),
    2731342: ("monodnaviria", "DNA Viruses"),
    2842242: ("ribozyviria", "RNA Viruses"),
    687329: ("anelloviridae", "DNA Viruses"),
    2559587: ("riboviria", "RNA Viruses"),
    9999999999: ("human viruses", "Viruses"),  # placeholder tax id
}  # excluded studies:
# "Johnson 2023",  # unpublished data
# "Cui 2023",  # untreated undigested sludge
# "Wang 2022",  # COVID-19 hospital wastewater
# "Petersen 2015",  # air plane waste
# "Hendriksen 2019",  # man hole"
# "Moritz 2019",  # university wastewater
# "Wu 2020",  # lung sample
# "Fierer 2022",  # university campus
# "Bohl 2022",


def assemble_data():
    plotting_data = []

    for study in studies:
        # Dropping studies that aren't WTP based
        if study not in [
            "Bengtsson-Palme 2016",
            "Munk 2022",
            "Brinch 2020",
            "Ng 2019",
            "Maritz 2019",
            "Brumfield 2022",
            "Rothman 2021",
            "Yang 2020",
            "Spurbeck 2023",
            "Crits-Christoph 2021",
        ]:
            continue

        for bioproject in metadata_papers[study]["projects"]:
            samples = metadata_bioprojects[bioproject]
            samples = metadata_bioprojects[bioproject]

            if study == "Bengtsson-Palme 2016":
                samples = [
                    sample
                    for sample in samples
                    if metadata_samples[sample]["fine_location"].startswith(
                        "Inlet"
                    )
                ]

            if study == "Ng 2019":
                samples = [
                    sample
                    for sample in samples
                    if metadata_samples[sample]["fine_location"] == "Influent"
                ]

            for sample in samples:
                modified_study = study
                if metadata_samples[sample].get("enrichment") == "panel":
                    continue
                if study == "Brumfield 2022":
                    print(metadata_samples[sample]["na_type"])
                    if metadata_samples[sample]["na_type"] == "DNA":
                        modified_study = "Brumfield 2022\n(DNA Subset)"
                        print("DNA subset was allocated")
                    else:
                        modified_study = "Brumfield 2022\n(RNA Subset)"
                        print("RNA subset was allocated")

                humanreads = "%s.humanviruses.tsv" % sample

                if not os.path.exists(f"../humanviruses/{humanreads}"):
                    subprocess.check_call(
                        [
                            "aws",
                            "s3",
                            "cp",
                            "s3://nao-mgs/%s/humanviruses/%s"
                            % (bioproject, humanreads),
                            "humanviruses/",
                        ]
                    )

                with open(f"../humanviruses/{humanreads}") as inf:
                    human_virus_reads = 0
                    human_virus_counts = {}
                    for line in inf:
                        (
                            line_taxid,
                            clade_assignments,
                            _,
                        ) = line.strip().split("\t")
                        clade_hits = int(clade_assignments)
                        line_taxid = int(line_taxid)
                        human_virus_reads += clade_hits
                        human_virus_counts[line_taxid] = clade_hits

                    human_virus_relative_abundance = (
                        human_virus_reads / metadata_samples[sample]["reads"]
                    )

                cladecounts = "%s.tsv.gz" % sample
                if not os.path.exists(f"../cladecounts/{cladecounts}"):
                    subprocess.check_call(
                        [
                            "aws",
                            "s3",
                            "cp",
                            "s3://nao-mgs/%s/cladecounts/%s"
                            % (bioproject, cladecounts),
                            "cladecounts/",
                        ]
                    )
                with gzip.open(f"../cladecounts/{cladecounts}") as inf:
                    taxa_abundances = {
                        "DNA Viruses": 0,
                        "RNA Viruses": 0,
                        "Viruses": 0,
                    }
                    for line in inf:
                        (
                            line_taxid,
                            _,
                            _,
                            clade_assignments,
                            _,
                        ) = line.strip().split()
                        taxid = int(line_taxid)
                        clade_hits = int(clade_assignments)
                        if taxid in target_taxa:
                            nucleic_acid_type = target_taxa[taxid][1]
                            relative_abundance = (
                                clade_hits / metadata_samples[sample]["reads"]
                            )

                            taxa_abundances[
                                nucleic_acid_type
                            ] += relative_abundance

                    plotting_data.append(
                        {
                            "study": modified_study,
                            "sample": sample,
                            "Human-Infecting Viruses": human_virus_relative_abundance,
                            **taxa_abundances,
                            **human_virus_counts,
                        }


    df = pd.Dataframe(plotting_data)
    

    return df


def load_taxonomic_data():
    parents = {}
    with open(os.path.join(dashboard, "nodes.dmp")) as inf:
        for line in inf:
            child_taxid, parent_taxid, child_rank, *_ = line.replace(
                "\t|\n", ""
            ).split("\t|\t")
            parent_taxid = int(parent_taxid)
            child_taxid = int(child_taxid)
            child_rank = child_rank.strip()
            parents[child_taxid] = (child_rank, parent_taxid)
    return parents


def get_family(taxid, parents):
    iteration_count = 0

    original_taxid = taxid

    current_taxid, current_rank, parent_taxid = (
        original_taxid,
    ) + parents.get(taxid, (None, None))

    while current_rank != "family" and parent_taxid is not None:
        current_taxid = parent_taxid
        current_rank, parent_taxid = parents.get(current_taxid, (None, None))
        iteration_count += 1

        if iteration_count > 100:
            break
    else:
        family_taxid = current_taxid
        return family_taxid


def get_taxid_name(target_taxid, taxonomic_names):
    tax_name = taxonomic_names[f"{target_taxid}"][0]
    return tax_name


def return_boxplot_plotting_df():
    df = assemble_box_plot_data()
    df_w_o_zeroes = df[(df.drop(["study", "sample"], axis=1) != 0).all(1)]
    df_plotting = df_w_o_zeroes.melt(
        id_vars=[
            "study",
            "sample",
        ],
        value_vars=[
            "Viruses",
            "RNA Viruses",
            "DNA Viruses",
            "Human-Infecting Viruses",
        ],
        var_name="Identifier",
        value_name="Relative Abundance",
    )

    return df_plotting


def return_barplot_plotting_df():
    df = assemble_bar_blot_data()

    parents = load_taxonomic_data()

    unique_numeric_cols = set(
        col for col in df.columns if isinstance(col, int)
    )

    species_family_mapping = {
        col: get_family(col, parents) for col in unique_numeric_cols
    }

    df.rename(columns=species_family_mapping, inplace=True)

    df_plotting = df.melt(
        id_vars=[
            "study",
            "sample",
        ],
        value_vars=[
            col
            for col in df.columns
            if col not in ["study", "sample"] and not pd.isna(col)
        ],
        var_name="family_taxid",
        value_name="read_count",
    )

    grouped_df = df_plotting.groupby(
        ["study", "family_taxid"]
    ).read_count.sum()

    print(grouped_df)

    grouped_df_w_o_zeroes = grouped_df[grouped_df != 0]

    df_pivot = grouped_df_w_o_zeroes.unstack(level=-1).fillna(0)

    df_normalized = df_pivot.div(df_pivot.sum(axis=1), axis=0)

    # taking the normal mean of three numbers:

    N_BIGGEST_FAMILIES = 9

    family_abundance_mean = df_normalized.apply(lambda x: np.mean(x), axis=0)

    top_families = family_abundance_mean.nlargest(N_BIGGEST_FAMILIES).index

    df_normalized["Other Viral Families"] = df_normalized.loc[
        :, ~df_normalized.columns.isin(top_families)
    ].sum(axis=1)

    df_normalized = df_normalized.loc[
        :, list(top_families) + ["Other Viral Families"]
    ]
    dict_family_name = {}

    for taxid in top_families.tolist():
        dict_family_name[taxid] = get_taxid_name(taxid, taxonomic_names)

    df_viral_family_names = df_normalized.rename(columns=dict_family_name)
    # return df_viral_family_names as csv
    df_viral_family_names.to_csv("viral_family_abundance.csv")

    return df_viral_family_names


def start():
    boxplot_df = return_boxplot_plotting_df()
    print(boxplot_df)
    barplot_df = return_barplot_plotting_df()
    print(barplot_df)

    study_nucleic_acid_mapping = {
        study: metadata["na_type"]
        for study, metadata in metadata_papers.items()
    }

    if "Brumfield 2022" in study_nucleic_acid_mapping:
        study_nucleic_acid_mapping["Brumfield 2022\n(DNA Subset)"] = "DNA"
        study_nucleic_acid_mapping["Brumfield 2022\n(RNA Subset)"] = "RNA"
        del study_nucleic_acid_mapping["Brumfield 2022"]

    print(study_nucleic_acid_mapping)

    boxplot_df = boxplot_df.reset_index()
    barplot_df = barplot_df.reset_index()

    boxplot_df["na_type"] = boxplot_df["study"].map(study_nucleic_acid_mapping)
    barplot_df["na_type"] = barplot_df["study"].map(study_nucleic_acid_mapping)

    na_type_order = ["DNA", "RNA"]

    boxplot_df = boxplot_df.sort_values(
        by="na_type",
        key=lambda col: col.map({k: i for i, k in enumerate(na_type_order)}),
    )

    barplot_df = barplot_df.sort_values(
        by="na_type",
        key=lambda col: col.map({k: i for i, k in enumerate(na_type_order)}),
    )

    boxplot_df["study"] = boxplot_df["study"].str.replace("-", "-\n")
    barplot_df["study"] = barplot_df["study"].str.replace("-", "-\n")

    fig = plt.figure(
        figsize=(9, 11),
    )

    gs = GridSpec(2, 2, height_ratios=[10, 5], figure=fig)

    ax0 = fig.add_subplot(gs[0, :])
    ax1 = fig.add_subplot(gs[1, :])

    sns.boxplot(
        data=boxplot_df,
        y="study",
        x="Relative Abundance",
        hue="Identifier",
        width=0.7,
        showfliers=False,
        ax=ax0,
    )

    y_order = ax0.get_yticklabels()

    print(y_order)

    y_order = [text.get_text() for text in y_order]

    ax0_title = ax0.set_title("a", fontweight="bold")
    ax0.set_xlabel("Relative Abundance among all reads")

    x, y = ax0_title.get_position()

    ax0_title.set_position((-0.18, 0))

    ax0.set_ylabel("")
    ax0.tick_params(left=False)

    ax0.set_xscale("log")
    ax0.set_xlim(right=1, left=1e-8)

    sns.despine(top=True, right=True, left=True, bottom=False)

    studies = boxplot_df["study"].unique()

    ax0.legend(
        loc=(0.0, -0.13),
        columnspacing=2.2,
        ncol=4,
        frameon=True,
        title="",
        fontsize=10,
    )

    for i in range(-7, 0):
        ax0.axvline(10**i, color="grey", linewidth=0.3, linestyle=":")

    for i in range(1, len(studies)):
        if i == 5:
            ax0.axhline(i - 0.5, color="black", linewidth=1, linestyle="-")

        else:
            ax0.axhline(i - 0.5, color="grey", linewidth=0.3, linestyle=":")

    ax0.text(1.15, -0.5, "DNA \nSequencing", ha="left")
    ax0.text(1.15, 4.7, "RNA \nSequencing", ha="left")

    ten_color_palette = [
        "#8dd3c7",
        "#f1c232",
        "#bebada",
        "#fb8072",
        "#80b1d3",
        "#fdb462",
        "#b3de69",
        "#fccde5",
        "#bc80bd",
        "#d9d9d9",
    ]

    barplot_df.set_index("study", inplace=True)

    barplot_df.loc[y_order].plot(
        kind="barh",
        stacked=True,
        color=ten_color_palette,
        ax=ax1,
    )

    # invert y axis
    ax1.invert_yaxis()

    ax1_title = ax1.set_title("b", fontweight="bold")

    x, y = ax1_title.get_position()

    ax1_title.set_position((-0.18, 0))

    ax1.set_xlabel("Relative Abundance among human-infecting viruses")

    ax1.tick_params(left=False)

    ax1.set_ylabel("")

    ax1.set_xlim(right=1, left=0)

    ax1.legend(
        loc=(0.0, -0.37),
        ncol=4,
        frameon=True,
        fontsize=9.6,
    )

    plt.tight_layout()
    plt.savefig("composite_fig_1.pdf", bbox_inches="tight")
    plt.savefig("composite_fig_1.png", dpi=600, bbox_inches="tight")


if __name__ == "__main__":
    start()
