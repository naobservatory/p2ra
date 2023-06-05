from pathogen_properties import *

background = """Hepatitis A is a vaccine-preventable liver-infection caused by
the hepatitis A virus. In the US, it's mostly spread by individual contact.
There is little seasonal variance in Hepatitis A incidence (https://www.cdc.
gov/vaccines/pubs/pinkbook/hepa.html). Viral shedding persists for 1 to 3
weeks."""


pathogen_chars = PathogenChars(
    na_type=NAType.RNA,
    enveloped=Enveloped.NON_ENVELOPED,
    # Using 12092 (Hepatitis A virus) instead of its child 208726 (Human
    # hepatitis A virus) because the MGS pipeline assigns reads to 12092.
    # Which happens because the Virus-Host DB
    # (https://www.genome.jp/virushostdb/view/) doesn't seem to know about
    # 208726.
    taxid=TaxID(12092),
)


us_incidence_absolute_2018 = IncidenceAbsolute(
    annual_infections=12474,
    confidence_interval=(17500, 27400),
    coverage_probability=0.95,
    country="United States",
    date="2018",
    source="https://www.cdc.gov/hepatitis/statistics/2018surveillance/HepA.htm",
)

us_population_2018 = Population(
    people=327.2 * 1e6,
    country="United States",
    date="2018",
    source="https://data.census.gov/table?q=2018+us+population&t=Civilian+Population",
)

acute_underreporting_factor = Scalar(
    scalar=2,
    confidence_interval=(1.4, 2.2),
    coverage_probability=0.95,
    country="United States",
    source="https://www.cdc.gov/hepatitis/statistics/2018surveillance/pdfs/2018HepSurveillanceRpt.pdf?#page=8",
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


def estimate_incidences() -> list[IncidenceRate]:
    us_2018 = us_incidence_absolute_2018.to_rate(us_population_2018)
    # 2019 was very similar to 2018, so extrapolation seems reasonable.
    us_2019 = dataclasses.replace(us_2018, date_source=Variable(date="2019"))
    # It's not clear how much the pandemic slowed down HepA transmission.
    # For now, assume somewhat dubiously that it didn't.
    us_2020 = dataclasses.replace(us_2018, date_source=Variable(date="2020"))
    us_2021 = dataclasses.replace(us_2018, date_source=Variable(date="2021"))

    return [
        us_2018,
        us_2019,
        us_2020,
        us_2021,
        king_county_confirmed_cases_rate_2017 * acute_underreporting_factor,
        king_county_confirmed_cases_rate_2018 * acute_underreporting_factor,
    ]


def estimate_prevalences() -> list[Prevalence]:
    return []
