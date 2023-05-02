from pathogen_properties import *

background = """According to the CDC, “Rhinoviruses are the most frequent 
cause of the common cold. In the United States children have an average of two 
rhinovirus infections each year, and adults have an average of one.”
Furthermore, Rhinovirus declined somewhat less than other respiratory viruses
during the pandemic (https://jamanetwork.com/journals/jamanetworkopen
fullarticle/2801017)."""


pathogen_chars = PathogenChars(
    na_type=NAType.RNA,
    enveloped=Enveloped.NON_ENVELOPED,
    taxid=TaxID(147711),
)

LA_county_child_proportion = Scalar(
    scalar=0.211,
    date=2021,
    source="https://www.census.gov/quickfacts/fact/table/losangelescountycalifornia#:~:text=Persons%20under-,18,-years%2C%20percent",
)

LA_county_adult_proportion = Scalar(
    scalar=1 - LA_county_child_proportion.scalar,
    date=2021,
    source="https://www.census.gov/quickfacts/fact/table/losangelescountycalifornia#:~:text=Persons%20under-,18,-years%2C%20percent",
)


def estimate_prevalences():
    return []
