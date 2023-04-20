from pathogen_properties import *

background = """Hepatitis A is a vaccine-preventable liver-infection caused by the hepatitis A virus. In the US,
it's mostly spread by individual contact. There is little seasonal variance in Hepatitis A incidence (https://www.cdc.gov/vaccines/pubs/pinkbook/hepa.html). Viral shedding persists for 1 to 3 weeks."""

pathogen_chars = PathogenChars(
    na_type=NAType.RNA,
    enveloped=Enveloped.NON_ENVELOPED,
    taxid=208726,
)

variables = {
    "us_incidence_absolute_2018": IncidenceAbsolute(
        annual_infections=12474,
        confidence_interval=(17500, 27400),
        country="United States",
        source="https://www.cdc.gov/hepatitis/statistics/2018surveillance/HepA.htm",
        start_date="2018-01-01",
        end_date="2018-12-31",
    ),
    "us_population_2018": Population(
        people=326.8 * 1e6,
        country="United States",
        source="https://www.cdc.gov/hepatitis/statistics/2018surveillance/HepA.htm",
        start_date="2018-01-01",
        end_date="2018-12-31",
    ),
    "hva_shedding_duration": SheddingDuration(
        days=14,
        confidence_interval=(7, 21),
        country="United States",
        source="https://www.cdc.gov/vaccines/pubs/pinkbook/hepa.html#:~:text=Viral%20shedding%20persists%20for%201%20to%203%20weeks.",
    ),
    "king_county_absolute_2017": IncidenceAbsolute(
        annual_infections=11,
        country="United States",
        state="Washington",
        county="King",
        start_date="2017-01-01",
        end_date="2017-12-31",
        source="https://doh.wa.gov/sites/default/files/2023-01/420-004-CDAnnualReport2021.pdf?uid=642c448518316#page=28",
    ),
    "king_county_absolute_2018": IncidenceAbsolute(
        annual_infections=14,
        country="United States",
        state="Washington",
        county="King",
        start_date="2018-01-01",
        end_date="2018-12-31",
        source="https://doh.wa.gov/sites/default/files/2023-01/420-004-CDAnnualReport2021.pdf?uid=642c448518316#page=28",
    ),
    "king_county_confirmed_cases_rate_2017": IncidenceRate(
        annual_infections_per_100k=0.5,
        country="United States",
        state="Washington",
        county="King",
        start_date="2017-01-01",
        end_date="2017-12-31",
        source="https://doh.wa.gov/sites/default/files/2023-01/420-004-CDAnnualReport2021.pdf?uid=642c448518316#page=28",
    ),
    "king_county_confirmed_cases_rate_2018": IncidenceRate(
        annual_infections_per_100k == 0.6,
        country="United States",
        state="Washington",
        county="King",
        start_date="2018-01-01",
        end_date="2018-12-31",
        source="https://doh.wa.gov/sites/default/files/2023-01/420-004-CDAnnualReport2021.pdf?uid=642c448518316#page=28",
    ),
}


def estimate_prevalences():
    return {
        "us_2018": variables["us_incidence_absolute_2018"]
        .to_rate(variables["us_population_2018"])
        .to_prevalence(variables["hva_shedding_duration"]),
        "king_county_2017": variables[
            "king_county_incidence_rate_2017"
        ].to_prevalence(variables["hva_shedding_duration"]),
        "king_county_2018": variables[
            "king_county_incidence_rate_2018"
        ].to_prevalence(variables["hva_shedding_duration"]),
    }