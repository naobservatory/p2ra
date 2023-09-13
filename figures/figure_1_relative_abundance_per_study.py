import os
import gzip
import json
import subprocess
import matplotlib.pyplot as plt  # type: ignore
import pandas as pd
import numpy as np
import seaborn as sns

os.chdir("../")


dashboard = os.path.expanduser("~/code/mgs-pipeline/dashboard/")

with open(os.path.join(dashboard, "human_virus_sample_counts.json")) as inf:
    human_virus_sample_counts = json.load(inf)

with open(os.path.join(dashboard, "metadata_samples.json")) as inf:
    metadata_samples = json.load(inf)

with open(os.path.join(dashboard, "metadata_bioprojects.json")) as inf:
    metadata_bioprojects = json.load(inf)

with open(os.path.join(dashboard, "metadata_papers.json")) as inf:
    metadata_papers = json.load(inf)

studies = list(metadata_papers.keys())

target_taxa = {
    10239: "viruses",
    2731341: "duplodnaviria",
    2559587: "riboviria",
    999999: "human viruses",  # placeholder taxid
}


def assemble_data():
    plotting_data = []
    for study in studies:
        if study == "Johnson 2023":
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

                humanreads = "%s.humanviruses.tsv" % sample

                if not os.path.exists(f"humanviruses/{humanreads}"):
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

                with open(f"humanviruses/{humanreads}") as inf:
                    human_virus_reads = 0
                    for line in inf:
                        (
                            line_taxid,
                            clade_assignments,
                            _,
                        ) = line.strip().split("\t")
                        clade_hits = int(clade_assignments)
                        human_virus_reads += int(clade_hits)
                    human_virus_relative_abundance = (
                        human_virus_reads / metadata_samples[sample]["reads"]
                    )

                cladecounts = "%s.tsv.gz" % sample
                if not os.path.exists(f"cladecounts/{cladecounts}"):
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
                with gzip.open(f"cladecounts/{cladecounts}") as inf:
                    taxa_abundances = {}
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
                        if taxid in list(target_taxa.keys()):
                            relative_abundance = (
                                clade_hits / metadata_samples[sample]["reads"]
                            )
                            taxa_abundances[
                                target_taxa[taxid]
                            ] = relative_abundance
                    plotting_data.append(
                        {
                            "study": study,
                            "sample": sample,
                            "human_viruses": human_virus_relative_abundance,
                            **taxa_abundances,
                        }
                    )
    df = pd.DataFrame(plotting_data)

    return df


def return_plotting_df():
    df = assemble_data()
    df_w_o_zeroes = df[(df.drop(["study", "sample"], axis=1) != 0).all(1)]
    df_plotting = df_w_o_zeroes.melt(
        id_vars=[
            "study",
            "sample",
        ],
        value_vars=["human_viruses", "viruses", "riboviria", "duplodnaviria"],
        var_name="Identifier",
        value_name="Relative Abundance",
    )
    study_nucleic_acid_mapping = {
        study: metadata["na_type"]
        for study, metadata in metadata_papers.items()
    }
    df_plotting["study"] = df_plotting["study"].astype("category")
    df_plotting["study"] = df_plotting["study"].cat.set_categories(
        sorted(
            df_plotting["study"].unique(),
            key=lambda x: study_nucleic_acid_mapping[x],
        ),
        ordered=True,
    )

    return df_plotting


def start():
    plotting_df = return_plotting_df()
    plt.figure(figsize=(8, 15))
    sns.boxplot(
        data=plotting_df,
        y="study",
        x="Relative Abundance",
        hue="Identifier",
        width=0.7,
        showfliers=False,
    )

    plt.ylabel("")
    plt.tick_params(left=False)

    plt.xscale("log")
    plt.xlim(right=1, left=1e-8)

    sns.despine(top=True, right=True, left=True, bottom=False)

    plt.legend(bbox_to_anchor=(1.1, 1), loc=2, borderaxespad=0.0)

    studies = plotting_df["study"].unique()
    for i in range(1, len(studies)):
        if i in [8, 11]:
            plt.axhline(i - 0.5, color="black", linewidth=1, linestyle="-")

        else:
            plt.axhline(i - 0.5, color="grey", linewidth=0.5, linestyle=":")

    plt.text(1.15, 7.7, "DNA+RNA", ha="left")
    plt.text(1.15, 10.7, "RNA", ha="left")
    plt.text(1.15, -0.5, "DNA", ha="left")

    plt.savefig(
        "figures/figure_1_relative_abundance_per_study.pdf",
        bbox_inches="tight",
    )
    print("hello")


if __name__ == "__main__":
    start()
