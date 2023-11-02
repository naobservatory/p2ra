import os
import gzip
import json
import subprocess
import matplotlib.pyplot as plt  # type: ignore
from matplotlib.gridspec import GridSpec
import matplotlib.ticker as ticker

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
}


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
                if metadata_samples[sample].get("enrichment") == "panel":
                    continue
                modified_study = study
                if study == "Brumfield 2022":
                    if metadata_samples[sample]["na_type"] == "DNA":
                        modified_study = "Brumfield 2022\n(DNA Subset)"
                    else:
                        modified_study = "Brumfield 2022\n(RNA Subset)"

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
                    human_virus_counts = {}
                    for line in inf:
                        (
                            line_taxid,
                            clade_assignments,
                            _,
                        ) = line.strip().split("\t")
                        clade_hits = int(clade_assignments)
                        line_taxid = int(line_taxid)
                        human_virus_counts[line_taxid] = clade_hits

                    plotting_data.append(
                        {
                            "study": modified_study,
                            "sample": sample,
                            **human_virus_counts,
                        }
                    )
    df = pd.DataFrame(plotting_data)

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


def return_formatted_data():
    df = assemble_data()

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
    df = return_formatted_data()

    caliciviridae_rothman = df.loc["Rothman 2021", "Caliciviridae"]
    caliciviridae_spurbeck = df.loc["Spurbeck 2023", "Caliciviridae"]
    polyomaviridae_munk = df.loc["Munk 2022", "Polyomaviridae"]
    polyomaviridae_brumfield = df.loc[
        "Brumfield 2022\n(DNA Subset)", "Polyomaviridae"
    ]

    print(
        f"the fraction of human-infecting virus reads mapping to Caliciviridae varied from {round(caliciviridae_rothman * 100)}% (Rothman 2022) to {round(caliciviridae_spurbeck * 100)}% (Spurbeck 2023), and Polyomaviridae relative abundance varied from {round(polyomaviridae_brumfield * 100)}% (Brumfield 2022 (DNA Subset)) to {round(polyomaviridae_munk*100)}% (Munk 2022)."
    )


if __name__ == "__main__":
    start()
