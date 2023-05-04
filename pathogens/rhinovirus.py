from pathogen_properties import *

background = """According to the CDC, “Rhinoviruses are the most frequent 
cause of the common cold. In the United States children have an average of two 
rhinovirus infections each year, and adults have an average of one.”
Furthermore, Rhinovirus declined somewhat less than other respiratory viruses
during the pandemic (https://jamanetwork.com/journals/jamanetworkopen
fullarticle/2801017).
Rhinovirus is not as seasonal as viruses like influenza or coronviruses.
(https://journals.plos.org/plosone/article/figure?id=10.1371/journal.pone.
0114871.g001)
(https://academic.oup.com/view-large/figure/89888454/195-6-773-fig003.jpeg)"""


pathogen_chars = PathogenChars(
    na_type=NAType.RNA,
    enveloped=Enveloped.NON_ENVELOPED,
    taxid=TaxID(147711, 147712, 463676),
)

LA_county_under_18_proportion = Scalar(
    scalar=0.211,
    date="2021",
    source="https://www.census.gov/quickfacts/fact/table/losangelescountycalifornia#:~:text=Persons%20under-,18,-years%2C%20percent",
)

LA_county_adult_proportion = Scalar(
    scalar=1 - LA_county_under_18_proportion.scalar,
    date="2021",
    source="https://www.census.gov/quickfacts/fact/table/losangelescountycalifornia#:~:text=Persons%20under-,18,-years%2C%20percent",
)

annual_rhinovirus_infections_children = Scalar(
    # rhinovirus_share_of_cold = 0.34
    scalar=3.7 * 0.34,
    source="https://sci-hub.ru/https://doi.org/10.1017/S0950268800050779#?page=6",
)

annual_rhinovirus_infections_adults = Scalar(
    # rhinovirus_share_of_cold = 0.34
    scalar=2 * 0.34,
    source="https://sci-hub.ru/https://doi.org/10.1017/S0950268800050779#?page=6",
)

annual_rhinovirus_infections_per_100k = IncidenceRate(
    annual_infections_per_100k=100_000
    * (
        LA_county_adult_proportion.scalar
        * annual_rhinovirus_infections_adults.scalar
        + LA_county_under_18_proportion.scalar
        * annual_rhinovirus_infections_children.scalar
    ),
)

rhinovirus_shedding_duration = SheddingDuration(
    days=12,
    confidence_interval=(10, 14),
    source="https://erj.ersjournals.com/content/44/1/169#:~:text=Virus%20shedding%20lasts%20on%20average%20for%2010%E2%80%9314%20days%20in%20immunocompetent%20subjects",
)

"""This article analyzes information from the following studies, and their 
analysis seems to make sense from an initial screening. https://www.rcgp.org.
uk/clinical-and-research/our-programmes/research-and-surveillance-centre/
public-health-data#:~:text=December%202022%20%C2%A0(ZIP%20file%2C%203.7%20MB
"""
pandemic_decrease_factor = Scalar(
    scalar=0.1,
    date="2020",
    source="https://www.bmj.com/content/370/bmj.m3182#:~:text=still%20around%20nine%20times%20fewer%20cases%20than%20the%20five%20year%20average%20for%20this%20time%20of%20year.",
)


def estimate_prevalences():
    return [
        annual_rhinovirus_infections_per_100k.to_prevalence(
            rhinovirus_shedding_duration
        ).scale(pandemic_decrease_factor)
    ]
