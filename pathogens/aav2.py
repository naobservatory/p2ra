import dataclasses

from pathogen_properties import *

background = """Adeno-associated virus 2 is a small replication-defective
virus. It belongs to the Dependoparvoviridae, a genus of viruses that is 
commonly used to construct viral vectors. AAV-2 has no clinical significance 
in humans, and seroprevalence studies are focused on patient groups that might 
receive a vector-based therapy. """

pathogen_chars = PathogenChars(
    na_type=NAType.DNA,
    enveloped=Enveloped.NON_ENVELOPED,
    taxid=TaxID(10804),
    selection=SelectionRound.ROUND_2,
)

seroprevalence_hemophilia_global_2021 = Prevalence(
    infections_per_100k=0.585 * 100_000,
    # Taking seropositivity from Figure 1A (un-weighted global seroprevalence)
    number_of_participants=513,
    # Though these participants are not representative of the general
    # population, hemophilia is not caused by AAV-2. Prevalence would thus be
    # at most affected by, e.g., lower socioeconomic status due to a higher
    # disease burden.
    country="Global",
    # Demographic composition:
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9063149/#:~:text=
    # (93K%2C%20docx)-,Supplemental%20data%3A,-Click%20here%20to
    date="2022",
    active=Active.LATENT,
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9063149/#:~:text=Seropositivity%20for%20(A)%20the%20global%20population"
    # This number matches AAV-2 seroprevalence in a study of 101 males with
    # Duchenne Muscular Dystrophy, showing a seroprevalence of 56%:
    # "https://pubmed.ncbi.nlm.nih.gov/36324212/#:~:text=We%20prospectively%20enrolled,and%20AAV8%20(47%25)."
)


def northern_european_average_seroprevalence() -> Prevalence:
    # Taking weighted average of Northern European AAV-5 seropositivity
    # numbers from Figure 1D, combined with participant numbers taken from the
    # supplement. We do this, because differences in seroprevalence between
    # countries are likely driven by small sample sizes, not by in-between-
    # country differences.
    participants_and_seroprevalence_by_country = {
        "France": (87, 0.605),
        "Germany": (90, 0.483),
        "United Kingdom": (17, 0.647),
    }
    prevalence_population_pairs: list[tuple[Prevalence, Population]] = []
    prevalence_population_pairs = [
        (
            Prevalence(
                infections_per_100k=seroprevalence * 100_000,
                number_of_participants=n_participants,
                country=country,
                date="2022",
                active=Active.LATENT,
                source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9063149/#:~:text=Seropositivity%20for%20(A)%20the%20global%20population",
            ),
            Population(
                people=n_participants,
                date="2022",
                country=country,
                source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9063149/#:~:text=(93K%2C%20docx)-,Supplemental%20data%3A,-Click%20here%20to",
            ),
        )
        for country, (
            n_participants,
            seroprevalence,
        ) in participants_and_seroprevalence_by_country.items()
    ]

    return Prevalence.weightedAverageByPopulation(*prevalence_population_pairs)


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

    # Extrapolating Northern European 2022 estimate to Denmark, backward in
    # time to 2015-2018:
    northern_europe_2022 = northern_european_average_seroprevalence()
    dk_2015 = dataclasses.replace(
        northern_europe_2022,
        date_source=Variable(date="2015"),
        location_source=Variable(country="Denmark"),
    )
    dk_2016 = dataclasses.replace(
        northern_europe_2022,
        date_source=Variable(date="2016"),
        location_source=Variable(country="Denmark"),
    )
    dk_2017 = dataclasses.replace(
        northern_europe_2022,
        date_source=Variable(date="2017"),
        location_source=Variable(country="Denmark"),
    )
    dk_2018 = dataclasses.replace(
        northern_europe_2022,
        date_source=Variable(date="2018"),
        location_source=Variable(country="Denmark"),
    )
    return [
        us_2020,
        us_2021,
        dk_2015,
        dk_2016,
        dk_2017,
        dk_2018,
    ]


def estimate_incidences():
    return []
