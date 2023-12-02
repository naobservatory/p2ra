import os
import json
import subprocess
import numpy as np
from scipy.stats import gmean  # type: ignore


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


def start():
    human_virus_shares = []
    zero_counts = 0

    for study in studies:
        # Dropping studies that aren't WTP based
        if study in [
            "Johnson 2023",  # unpublished data
            "Cui 2023",  # untreated undigested sludge
            "Wang 2022",  # COVID-19 hospital wastewater
            "Petersen 2015",  # air plane waste
            "Hendriksen 2019",  # man hole"
            "Moritz 2019",  # university wastewater
            "Wu 2020",  # lung sample
            "Fierer 2022",  # university campus
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

                    # print (human_virus_relative_abundance)

                    if human_virus_relative_abundance == 0.0:
                        zero_counts += 1
                        continue

                    human_virus_shares.append(human_virus_relative_abundance)
    # Dropping all zeros from human_virus_shares

    perc_zero_human_read_samples = round(
        (zero_counts / len(human_virus_shares)) * 100, 2
    )

    gmean_viral_share = round(gmean(human_virus_shares), 7)

    return f"When dropping {zero_counts} samples without human reads ({perc_zero_human_read_samples}% of all samples), the geometric mean of samples' human read share is {gmean_viral_share * 100}% \n Put differently 1 in {round(1 / gmean_viral_share)} reads in the wastewater samples are human reads"


if __name__ == "__main__":
    print(start())
