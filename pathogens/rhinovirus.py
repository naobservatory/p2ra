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
    date="2021",
    source="https://www.census.gov/quickfacts/fact/table/losangelescountycalifornia#:~:text=Persons%20under-,18,-years%2C%20percent",
)

LA_county_adult_proportion = Scalar(
    scalar=1 - LA_county_child_proportion.scalar,
    date="2021",
    source="https://www.census.gov/quickfacts/fact/table/losangelescountycalifornia#:~:text=Persons%20under-,18,-years%2C%20percent",
)

annual_rhinovirus_infections_children = Scalar(
    scalar=2,
    source="https://www.cdc.gov/ncird/rhinoviruses-common-cold.html",
)

annual_rhinovirus_infections_adults = Scalar(
    scalar=1,
    source="https://www.cdc.gov/ncird/rhinoviruses-common-cold.html",
)

annual_rhinovirus_infections_per_100k = IncidenceRate(
    annual_infections_per_100k=100000
    * (
        LA_county_adult_proportion.scalar
        * annual_rhinovirus_infections_adults.scalar
        + LA_county_child_proportion.scalar
        * annual_rhinovirus_infections_children.scalar
    ),
)

rhinovirus_shedding_duration = SheddingDuration(
    days=8.5,
    source="https://erj.ersjournals.com/content/44/1/169#:~:text=Virus%20shedding%20lasts%20on%20average%20for%2010%E2%80%9314%20days%20in%20immunocompetent%20subjects",
)


def estimate_prevalences():
    return [
        annual_rhinovirus_infections_per_100k.to_prevalence(
            rhinovirus_shedding_duration
        )
    ]
