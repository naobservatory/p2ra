import csv

from pathogen_properties import *
from populations import us_population

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

cdc_underreporting_factor_2019 = Scalar(
    scalar=2,
    confidence_interval=(1.4, 2.2),
    coverage_probability=0.95,
    country="United States",
    source="https://www.cdc.gov/hepatitis/statistics/2019surveillance/Introduction.htm#Technical:~:text=The%20published%20multipliers%20have%20since%20been%20corrected%20by%20CDC%20to%20indicate%20that%20each%20reported%20case%20of%20acute%20hepatitis%20A%20represents%202.0%20estimated%20infections%20(95%25%20bootstrap%20CI%3A%201.4%E2%80%932.2)",
)

us_population_2018 = Population(
    people=327.2 * 1e6,
    country="United States",
    date="2018",
    source="https://data.census.gov/table?q=2018+us+population&t=Civilian+Population",
)

us_population_2019 = Population(
    people=328.2 * 1e6,
    country="United States",
    date="2019",
    source="https://data.census.gov/table?q=2019+us+population&t=Civilian+Population&tid=ACSDP1Y2019.DP05",
)

us_estimated_incidence_absolute_2018 = IncidenceAbsolute(
    annual_infections=12_474,
    confidence_interval=(17_500, 27_400),
    coverage_probability=0.95,
    country="United States",
    date="2018",
    source="https://www.cdc.gov/hepatitis/statistics/2018surveillance/HepA.htm",
)

us_estimated_incidence_absolute_2019 = IncidenceAbsolute(
    annual_infections=37_700,
    confidence_interval=(26_400, 41_500),
    coverage_probability=0.95,
    country="United States",
    date="2019",
    source="https://www.cdc.gov/hepatitis/statistics/2019surveillance/Introduction.htm#Technical:~:text=During%202019%2C%20a%20total%20of%2018%2C846%20hepatitis%20A%20cases%20were%20reported%20to%20CDC%2C%20corresponding%20to%2037%2C700%20estimated%20infections%20(95%25%20confidence%20interval%20%5BCI%5D%3A%2026%2C400%E2%80%9341%2C500)%20after%20adjusting%20for%20case%20underascertainment%20and%20underreporting%20(see%20Technical%20Notes)%20(9).",
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

ohio_reported_incidence_rate_2020 = IncidenceRate(
    annual_infections_per_100k=2.4,
    date="2020",
    country="United States",
    state="Ohio",
    source="https://www.cdc.gov/hepatitis/statistics/2020surveillance/hepatitis-a/table-1.1.htm#:~:text=Ohio,2.4",
    # Source provides yearly statewide data going back to 2016 for 30+ states
)

ohio_county_hav_incidences = {}


# Data source: https://odh.ohio.gov/know-our-programs/Hepatitis-Surveillance-Program/Hepatitis-A-Statewide-Community-Outbreak/
# Data is incidence over 4.5 year period
with open(prevalence_data_filename("havCaseCountsOhioCounties.csv")) as file:
    reader = csv.reader(file)
    next(reader)
    for row in reader:
        row[0] = row[0] + " County"
        ohio_county_hav_incidences[row[0]] = IncidenceAbsolute(
            country="United States",
            state="Ohio",
            date="2021",
            county=row[0],
            annual_infections=int(row[1]) / 4.5,
        ).to_rate(us_population(year=2021, state="Ohio", county=row[0]))


def estimate_incidences() -> list[IncidenceRate]:
    estimates = [
        us_estimated_incidence_absolute_2018.to_rate(us_population_2018),
        us_estimated_incidence_absolute_2019.to_rate(us_population_2019),
        king_county_confirmed_cases_rate_2017 * cdc_underreporting_factor_2019,
        king_county_confirmed_cases_rate_2018 * cdc_underreporting_factor_2019,
        ohio_reported_incidence_rate_2020,
    ]
    for item in ohio_county_hav_incidences:
        estimates.append((ohio_county_hav_incidences[item]))
    return estimates


def estimate_prevalences() -> list[Prevalence]:
    return []
