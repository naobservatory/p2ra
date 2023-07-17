#!/usr/bin/env python3
import pathogen_properties
import pathogens


def start():
    taxid_to_name = {}
    for pathogen_name, pathogen in pathogens.pathogens.items():
        for taxids in pathogen_properties.by_taxids(
            pathogen.pathogen_chars, pathogen.estimate_prevalences()
        ):
            for taxid in taxids:
                taxid_to_name[taxid] = pathogen_name
        for taxids in pathogen_properties.by_taxids(
            pathogen.pathogen_chars, pathogen.estimate_incidences()
        ):
            for taxid in taxids:
                taxid_to_name[taxid] = pathogen_name

    print("taxid", "filename", "human_readable", sep="\t")
    for taxid, name in sorted(taxid_to_name.items()):
        print(taxid, name, pathogens.tidy_name(name, [taxid]), sep="\t")


if __name__ == "__main__":
    start()
