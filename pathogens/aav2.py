import dataclasses

from pathogen_properties import *

background = """Adeno-associated virus 2 is a small replication-defective
virus. It belongs to the Dependoparvovirus, a genus of viruses that is commonly used to construct viral vectors. AAV-2 has no clinical significance in humans, and seroprevalence studies are focused on patient groups that might receive a vector-based therapy. """

pathogen_chars = PathogenChars(
    na_type=NAType.DNA,
    enveloped=Enveloped.NON_ENVELOPED,
    taxid=TaxID(10804),
)

seroprevalence_hemophilia_global_2021 = Prevalence(
    infections_per_100k=0.585 * 100_000,
    # "Factoring in the prevalence of HA in the countries being assayed, the
    # global weighted average of AAV5 seroprevalence in people with HA was 29.
    # 7%. For other AAV serotypes, global seroprevalence was 58.5% for AAV2,
    # [...]"
    number_of_participants=513,
    # Though these participants are not representative of the general
    # population, hemophilia is not caused by AAV-2. Prevalence would thus be
    # at most affected by, e.g., lower socioeconomic status due to a higher
    # disease burden.
    country="Global",
    # Demographic composition:
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9063149/#:~:text=
    # (93K%2C%20docx)-,Supplemental%20data%3A,-Click%20here%20to
    date="2021",
    active=Active.LATENT,
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9063149/#:~:text=Factoring%20in%20the,44.7%25%20for%20AAVrh10."
    # This number matches AAV-2 seroprevalence in a study of 101 males with
    # Duchenne Muscular Dystrophy, showing a seroprevalence of 56%:
    # "https://pubmed.ncbi.nlm.nih.gov/36324212/#:~:text=We%20prospectively%20enrolled,and%20AAV8%20(47%25)."
)


def estimate_prevalences() -> list[Prevalence]:
    # We assume that global seroprevalence will be similar to seroprevalence
    # in the US. This is also what we find in the US-American participants of
    # the same study (n=71, prevalence(AAV-2)=53.5%, found in Figure 1C.
    # Source:
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9063149/#:~:text=for%20(C)-,AAV2,-%2C%20(D)
    #  We also assume that seroprevalence remains constant over time, given
    # that AAV-2 is not known to cause disease, and thus isn't being treated
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
