#!/usr/bin/env python3

import mgs
from pathogens import pathogens

if __name__ == "__main__":
    repo = mgs.GitHubRepo(
        user="naobservatory",
        repo="mgs-pipeline",
        branch="47e2025f35168d3f414ae62928f6a14dd3f7c23d",
    )
    bp_data = mgs.load_bioprojects(repo)
    sample_data = mgs.load_sample_attributes(repo)
    counts = mgs.load_sample_counts(repo)
    taxtree = mgs.load_tax_tree(repo)

    bioproject = mgs.BioProject("PRJNA729801")  # Rothman
    samples = bp_data[bioproject]
    sample_attribs = {s: sample_data[s] for s in samples}
    fine_locs = set(sample_attribs[s].fine_location for s in samples)

    for pathogen_name, pathogen in pathogens.items():
        print(pathogen_name)
        taxid = pathogen.pathogen_chars.taxid
        subtree = taxtree[taxid]
        if not subtree:
            continue
        virus_counts = mgs.count_reads(subtree, counts)
        print(" All", sum(virus_counts[s] for s in samples), sep="\t")
        for fine_loc in fine_locs:
            print(
                f" {fine_loc}",
                sum(
                    virus_counts[s]
                    for s in samples
                    if sample_attribs[s].fine_location == fine_loc
                ),
                sep="\t",
            )
