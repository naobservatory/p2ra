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

incidence_underreporting_scalar = Scalar(
    scalar=1 / 0.59,
    confidence_interval=(
        1 / 0.84,
        1 / 0.32,
    ),
    coverage_probability=0.95,
    country="United States",
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4906888/#:~:text=Diagnosed%20hepatitis%20A,clear%20reporting%20responsibilities.",
)

king_county_absolute_2017 = IncidenceAbsolute(
    annual_infections=11,
    country="United States",
    state="Washington",
    county="King County",
    date="2017",
    source="https://doh.wa.gov/sites/default/files/2023-01/420-004-CDAnnualReport2021.pdf?uid=642c448518316#page=28",
)

king_county_absolute_2018 = IncidenceAbsolute(
    annual_infections=14,
    country="United States",
    state="Washington",
    county="King County",
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

ohio_hav_prevalence_2020 = Prevalence(
    infections_per_100k=2.4,
    date="2020",
    country="United States",
    state="Ohio",
    active=Active.ACTIVE,
    source="https://www.cdc.gov/hepatitis/statistics/2020surveillance/hepatitis-a/figure-1.3.htm",
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
        us_incidence_absolute_2018.to_rate(us_population_2018),
        king_county_confirmed_cases_rate_2017
        * incidence_underreporting_scalar,
        king_county_confirmed_cases_rate_2018
        * incidence_underreporting_scalar,
    ]
    for item in ohio_county_hav_incidences:
        estimates.append((ohio_county_hav_incidences[item]))
    return estimates


def estimate_prevalences() -> list[Prevalence]:
    return [ohio_hav_prevalence_2020]
