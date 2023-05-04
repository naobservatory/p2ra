#!/usr/bin/env python3

import mgs
from pathogens import pathogens

if __name__ == "__main__":
    bioproject = mgs.BioProject("PRJNA729801")  # Rothman
    mgs_data = mgs.MGSData.from_repo()
    samples = mgs_data.sample_attributes(
        bioproject, enrichment=mgs.Enrichment.VIRAL
    )

    fine_locs = set(attribs.fine_location for _, attribs in samples.items())

    for pathogen_name, pathogen in pathogens.items():
        print(pathogen_name)
        taxid = pathogen.pathogen_chars.taxid
        virus_reads = mgs_data.viral_reads(bioproject, taxid)
        print(" All", sum(virus_reads[s] for s in samples), sep="\t")
        for fine_loc in fine_locs:
            print(
                f" {fine_loc}",
                sum(
                    virus_reads[s]
                    for s, attribs in samples.items()
                    if attribs.fine_location == fine_loc
                ),
                sep="\t",
            )
