from p2ra_prevalences.pathogen_properties import *


background = """Hepatitis A is a vaccine-preventable liver-infection caused by the hepatitis A virus. In the US,
it's mostly spread by individual contact. There is little seasonal variance in Hepatitis A incidence (https://www.cdc.gov/vaccines/pubs/pinkbook/hepa.html). Viral shedding persists for 1 to 3 weeks."""

pathogen_chars = PathogenChars(
    na_type=NAType.RNA,
    enveloped=Enveloped.NON_ENVELOPED,
    taxid=208726,
)

prevalence_est = {  # prevalence estimators
    "us_incidence_absolute_2018": PrevalenceEstimator(
        value=12474,
        unit="cases",
        value_type="incidence_absolute",
        country="United States",
        source="https://www.cdc.gov/hepatitis/statistics/2018surveillance/index.htm",
        start_date="2018-01-01",
        end_date="2018-12-31",
    ),
    "us_population_2018": PrevalenceEstimator(
        value=326.8 * 1e6,
        unit="people",
        value_type="population",
        country="United States",
        source="https://www.google.com/search?q=us+population+2018",
        start_date="2018-01-01",
        end_date="2018-12-31",
    ),
    "hva_shedding_duration": PrevalenceEstimator(
        value=14,
        confidence_interval=(7, 21),
        unit="days",
        value_type="shedding_duration",
        country="United States",
        source="https://www.cdc.gov/vaccines/pubs/pinkbook/hepa.html#:~:text=Viral%20shedding%20persists%20for%201%20to%203%20weeks.",
    ),
    "king_county_confirmed_cases_rate_2017": PrevalenceEstimator(
        value=0.5,
        unit="cases_per_100k",
        value_type="incidence_rate",
        country="United States",
        state="Washington",
        county="King",
        start_date="2017-01-01",
        end_date="2017-12-31",
        source="https://doh.wa.gov/sites/default/files/2023-01/420-004-CDAnnualReport2021.pdf?uid=642c448518316#page=28",
    ),
    "king_county_confirmed_cases_rate_2018": PrevalenceEstimator(
        value=0.6,
        unit="cases_per_100k",
        value_type="incidence_rate",
        country="United States",
        state="Washington",
        county="King",
        start_date="2018-01-01",
        end_date="2018-12-31",
        source="https://doh.wa.gov/sites/default/files/2023-01/420-004-CDAnnualReport2021.pdf?uid=642c448518316#page=28",
    ),
}


def calculate_prevalence_from_incidence_rate(
    incidence_rate: float, shedding_duration: int
) -> float:
    prevalence = incidence_rate * (shedding_duration / 365)
    return prevalence


def calculate_prevalence_from_absolute_incidence(
    absolute_incidence: int, population: int, shedding_duration: int
) -> float:

    incidence_rate = (absolute_incidence / population) * 100000

    prevalence = incidence_rate * (shedding_duration / 365)

    return prevalence


prevalence_vars = {
    "us_prevalence_2018": calculate_prevalence_from_absolute_incidence(
        prevalence_est["us_incidence_absolute_2018"].value,
        prevalence_est["us_population_2018"].value,
        prevalence_est["hva_shedding_duration"].value,
    ),
    "king_county_prevalence_2017": calculate_prevalence_from_incidence_rate(
        prevalence_est["king_county_incidence_rate_2017"].value,
        prevalence_est["hva_shedding_duration"].value,
    ),
    "king_county_prevalence_2018": calculate_prevalence_from_incidence_rate(
        prevalence_est["king_county_incidence_rate_2018"].value,
        prevalence_est["hva_shedding_duration"].value,
    ),
}
