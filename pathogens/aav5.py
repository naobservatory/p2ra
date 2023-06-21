import dataclasses

from pathogen_properties import *

background = """Adeno-associated virus 5 is a small replication-defective
virus. It belongs to the Dependoparvoviridae, a genus of viruses that is 
commonly used to construct viral vectors. AAV-5 has no known clinical 
significance in humans, and seroprevalence studies are focused on patient 
groups that might receive a vector-based therapy. """

pathogen_chars = PathogenChars(
    na_type=NAType.DNA,
    enveloped=Enveloped.NON_ENVELOPED,
    taxid=TaxID(82300),
)

seroprevalence_hemophilia_global_2021 = Prevalence(
    infections_per_100k=0.446 * 100_000,
    # "Factoring in the prevalence of HA in the countries being assayed, the
    # global weighted average of AAV5 seroprevalence in people with HA was 29.
    # 7%."
    number_of_participants=513,
    # Though these participants are not representative of the general
    # population, hemophilia is not caused by AAV-5. Prevalence would thus be
    # at most affected by, e.g., lower socioeconomic status due to a higher
    # disease burden.
    country="Global",
    # Demographic composition:
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9063149/#:~:text=
    # (93K%2C%20docx)-,Supplemental%20data%3A,-Click%20here%20to
    date="2021",
    active=Active.LATENT,
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9063149/#:~:text=Factoring%20in%20the,44.7%25%20for%20AAVrh10."
    # As observed in AAV2 and AAV6, seroprevalence numbers of this study 
    # broadly agree with those in other studies


def estimate_prevalences() -> list[Prevalence]:
    # We assume that global seroprevalence will be similar to seroprevalence
    # in the US. This is also what we find in the US-American participants of
    # the same study (n=71, prevalence(AAV-6)=38.0%. Source:
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9063149/#:~:text=AAV5%2C%20(E)-,AAV6,-%2C%20(F)
    #  We also assume that seroprevalence remains constant over time, given
    # that AAV-6 is not known to cause disease, and thus isn't being treated
    # or vaccinated against.
    us_2020 = dataclasses.replace(
        seroprevalence_hemophilia_global_2021,
        date_source=Variable(date="2020"),
        location_source=Variable(country="United States"),
    )
    us_2021 = dataclasses.replace(
        us_2020,
        date_source=Variable(date="2021"),
    )
    return [
        us_2020,
        us_2021,
    ]


def estimate_incidences():
    return []
