from pathogen_properties import *

background = """Hepatitis A is a vaccine-preventable liver-infection caused by
the hepatitis A virus. In the US, it's mostly spread by individual contact.
There is little seasonal variance in Hepatitis A incidence (https://www.cdc.
gov/vaccines/pubs/pinkbook/hepa.html). Viral shedding persists for 1 to 3
weeks."""


pathogen_chars = PathogenChars(
    na_type=NAType.RNA,
    enveloped=Enveloped.NON_ENVELOPED,
    taxid=208726,
)


us_incidence_absolute_2018 = IncidenceAbsolute(
    annual_infections=12474,
    confidence_interval=(17500, 27400),
    coverage_probability_perc_point=95,
    country="United States",
    source="https://www.cdc.gov/hepatitis/statistics/2018surveillance/HepA.htm",
    date="2018",
)

us_population_2018 = Population(
    people=327.2 * 1e6,
    country="United States",
    source="https://data.census.gov/table?q=2018+us+population&t=Civilian+Population",
    date="2018",
)

hva_shedding_duration = SheddingDuration(
    days=14,
    confidence_interval=(7, 21),
    country="United States",
    source="https://www.cdc.gov/vaccines/pubs/pinkbook/hepa.html#:~:text=Viral%20shedding%20persists%20for%201%20to%203%20weeks.",
)

incidence_underreporting_scalar = Scalar(
    scalar=1 / 0.59,
    confidence_interval=(
        1 / 0.84,
        1 / 0.32,
    ),
    coverage_probability_perc_point=95,
    country="United States",
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4906888/#:~:text=Diagnosed%20hepatitis%20A,clear%20reporting%20responsibilities.",
)

king_county_absolute_2017 = IncidenceAbsolute(
    annual_infections=11,
    country="United States",
    state="Washington",
    county="King",
    date="2017",
    source="https://doh.wa.gov/sites/default/files/2023-01/420-004-CDAnnualReport2021.pdf?uid=642c448518316#page=28",
)

king_county_absolute_2018 = IncidenceAbsolute(
    annual_infections=14,
    country="United States",
    state="Washington",
    county="King",
    date="2018",
    source="https://doh.wa.gov/sites/default/files/2023-01/420-004-CDAnnualReport2021.pdf?uid=642c448518316#page=28",
)

king_county_confirmed_cases_rate_2017 = IncidenceRate(
    annual_infections_per_100k=0.5,
    country="United States",
    state="Washington",
    county="King",
    date="2017",
    source="https://doh.wa.gov/sites/default/files/2023-01/420-004-CDAnnualReport2021.pdf?uid=642c448518316#page=28",
)


king_county_confirmed_cases_rate_2018 = IncidenceRate(
    annual_infections_per_100k=0.6,
    country="United States",
    state="Washington",
    county="King",
    date="2018",
    source="https://doh.wa.gov/sites/default/files/2023-01/420-004-CDAnnualReport2021.pdf?uid=642c448518316#page=28",
)


def estimate_prevalences():
    return [
        us_incidence_absolute_2018.to_rate(us_population_2018).to_prevalence(
            hva_shedding_duration
        ),
        king_county_confirmed_cases_rate_2017.to_prevalence(
            hva_shedding_duration
        ).scale(incidence_underreporting_scalar),
        king_county_confirmed_cases_rate_2018.to_prevalence(
            hva_shedding_duration
        ).scale(incidence_underreporting_scalar),
    ]
