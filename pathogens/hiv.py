import dataclasses

from pathogen_properties import *
from populations import us_population

background = """HIV is a sexually-transmitted retrovirus which gradually
weakens the immune system."""


pathogen_chars = PathogenChars(
    na_type=NAType.RNA,
    enveloped=Enveloped.ENVELOPED,
    taxid=TaxID(11676),
)

# Controlled HIV has basically no sheddding, and we're really only interested
# in uncontrolled HIV.  The CDC estimates that 57% of Americans with HIV were
# virally suppressed in 2019.
us_unsuppressed_fraction_2019 = Scalar(
    scalar=1 - 0.57,
    country="United States",
    date="2019",
    source="https://www.cdc.gov/hiv/library/reports/hiv-surveillance/vol-26-no-2/content/national-profile.html#:~:text=57%25%20were%20virally%20suppressed",
)

us_infected_2019 = PrevalenceAbsolute(
    infections=1.2e6,
    country="United States",
    date="2019",
    active=Active.LATENT,
    source="https://www.cdc.gov/hiv/library/reports/hiv-surveillance/vol-26-no-2/content/national-profile.html#:~:text=Among%20the%20estimated-,1.2%20million%20people,-living%20with%20HIV",
)

us_population_2019 = Population(
    people=328_231_337,
    country="United States",
    date="2019-01-01",
    source="https://www.census.gov/newsroom/press-releases/2019/new-years-population.html",
)

la_unsuppressed_fraction_2020 = Scalar(
    scalar=1 - 0.6,
    country="United States",
    state="California",
    county="Los Angeles County",
    date="2020",
    source="https://web.archive.org/web/20201202004910/https://www.lacounty.hiv/",
)

la_infected_2020 = PrevalenceAbsolute(
    infections=57_700,
    country="United States",
    state="California",
    county="Los Angeles County",
    date="2020",
    active=Active.LATENT,
    source="https://web.archive.org/web/20201202004910/https://www.lacounty.hiv/",
)


def estimate_prevalences() -> list[Prevalence]:
    us_2019 = us_infected_2019.to_rate(us_population_2019) * us_unsuppressed_fraction_2019

    # HIV prevalence should be close to constant, so it's fine to figure that
    # it's about the same in 2020 and 2021 as it was in 2019.
    us_2020 = dataclasses.replace(
        us_2019, date_source=Variable(date="2020"))
    us_2021 = dataclasses.replace(
        us_2019, date_source=Variable(date="2021"))
    return [
        us_2019,
        us_2020,
        us_2021,
        la_infected_2020.to_rate(
            us_population(
                state="California", county="Los Angeles County", year=2020
            )
        )
        * la_unsuppressed_fraction_2020,
    ]


def estimate_incidences() -> list[IncidenceRate]:
    return []
