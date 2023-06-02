import csv
from collections import Counter

from pathogen_properties import *
from populations import us_population

background = """Epstein-Barr virus (EBV) is a type of  herpes virus that
infects humans and is known to cause infectious mononucleosis, also known
as mono or glandular fever. EBV is a common virus that is transmitted 
through contact with infected saliva, such as through kissing, sharing 
utensils, or close contact with an infected person's respiratory droplets.
EBV is a widespread virus that can persist in the body for life, although
most people infected with EBV do not develop symptoms or have mild symptoms
that resemble the flu."""

# TODO: Simon will look into incidence of mononucleosis in the US. There is
# also the question of mono being a proxy for EBV infection, or activation.

pathogen_chars = PathogenChars(
    na_type=NAType.DNA,
    enveloped=Enveloped.ENVELOPED,
    taxid=TaxID(10376),
)

us_fraction_u18 = Scalar(
    scalar=0.222,
    source="https://www.census.gov/quickfacts/fact/table/US/LFE046221#:~:text=%EE%A0%BF-,22.2%25,-Persons%2065%20years",
)
us_population_u18 = us_population(year=2022) * us_fraction_u18

us_population_18plus = us_population(year=2022) - us_population_u18

under_18_population_US = Population(
    people=333_287_557 * 0.222,
    # This is the same number as is used in populations.py
    date="2022",
    country="United States",
    source="",
)

over_18_population_US = Population(
    people=333_287_557 * (1 - 0.222),
    country="United States",
    date="2022",
    source="https://www.census.gov/quickfacts/fact/table/US/LFE046221#:~:text=%EE%A0%BF-,22.2%25,-Persons%2065%20years",
)

uk_seroprevalence_0_to_25 = Prevalence(
    infections_per_100k=0.853 * 100_000,
    # Given very high seroprevalence, and lifetime persistence of EBV, we
    # we treat this data as corresponding to seroprevalence across the
    # entire population.
    country="United Kingdom",
    start_date="2002",
    end_date="2013",
    number_of_participants=2325,
    active=Active.LATENT,
    source="https://bmcpublichealth.biomedcentral.com/articles/10.1186/s12889-020-09049-x#:~:text=1982/2325%20individuals%20(85.3%25)%20were%20EBV%20seropositive",
)

nhanes_6_19_yo_seroprevalence_estimate_2003_2010 = Prevalence(
    infections_per_100k=0.665 * 100_000,
    confidence_interval=(0.643 * 100_000, 0.687 * 100_000),
    coverage_probability=0.95,
    country="United States",
    start_date="2003",
    end_date="2010",
    number_of_participants=8417,
    active=Active.LATENT,
    source="https://pubmed.ncbi.nlm.nih.gov/23717674/#:~:text=Overall%20EBV%20seroprevalence%20was%2066.5%25%20(95%25%20CI%2064.3%25%2D68.7%25.)",
)

nhanes_18_19_yo_seroprevalence_estimate_2009_2010 = Prevalence(
    # Only using data from NHANES 2009-2010, for 18-19 year olds.
    infections_per_100k=0.89 * 100_000,
    # Given very high seroprevalence, and lifetime persistence of EBV, we
    # we treat this data as corresponding to seroprevalence across the
    # entire adult population.
    confidence_interval=(0.81 * 100_000, 0.94 * 100_000),
    coverage_probability=0.95,
    number_of_participants=508,
    start_date="2009",
    end_date="2010",
    country="United States",
    active=Active.LATENT,
    source="https://academic.oup.com/jid/article/208/8/1286/2192838#:~:text=273-,89%20(81%E2%80%9394)D%C2%A0,-.337%C2%A0",
)


us_seroprevalence_2003_2010 = Prevalence.weightedAverageByPopulation(
    (nhanes_6_19_yo_seroprevalence_estimate_2003_2010, under_18_population_US),
    (nhanes_18_19_yo_seroprevalence_estimate_2009_2010, over_18_population_US),
    # As noted previously, given very high seroprevalence, and lifetime
    # persistence of EBV, we treat this 18 to 19yo data as corresponding
    # to seroprevalence across the entire adult population.
)


# Source for CSV: https://www.census.gov/data-tools/demo/idb/#/pop?COUNTRY_YEAR=2023&COUNTRY_YR_ANIM=2023&FIPS_SINGLE=DA&menu=popViz&FIPS=DA&POP_YEARS=2023&popPages=BYAGE
with open(prevalence_data_filename("DenmarkPopulationData.csv")) as file:
    denmark_age_groups: Counter[tuple[int, int]] = Counter()
    # (min_age, max_age) -> people within that age range
    reader = csv.reader(file)
    next(reader)  # Skip the header

    for row in reader:
        # Skip the total population row
        if row[4] == "TOTAL":
            continue

        age = int(
            row[4].replace("+", "")
        )  # Extract the age from the first column & turn "100+" into 100

        # Remove spaces from the numbers and add them to the respective
        # age group
        population = int(row[5].replace(" ", ""))

        if age == 0:
            denmark_age_groups[(0, 0)] += population
        elif 1 <= age <= 3:
            denmark_age_groups[(1, 3)] += population
        elif 4 <= age <= 14:
            denmark_age_groups[(4, 14)] += population
        elif 15 <= age <= 17:
            denmark_age_groups[(15, 17)] += population
        elif 18 <= age <= 29:
            denmark_age_groups[(18, 29)] += population
        elif 30 <= age <= 100:
            denmark_age_groups[(30, 100)] += population
        else:
            assert False

    denmark_populations = {}
    for (min_age, max_age), people in denmark_age_groups.items():
        denmark_populations[min_age, max_age] = Population(
            people=people,
            date="2023",
            country="Denmark",
            source="https://www.census.gov/data-tools/demo/idb/#/pop?COUNTRY_YEAR=2023&COUNTRY_YR_ANIM=2023&FIPS_SINGLE=DA&menu=popViz&FIPS=DA&POP_YEARS=2023&popPages=BYAGE",
        )

DENMARK_SEROPREVALENCE_SOURCE = "10.3109/inf.1983.15.issue-4.03"

denmark_seroprevalences = {}
for age_range, seroprevalence in [
    ((0, 0), 0.15 * 100_000),
    ((1, 3), 0.28 * 100_000),
    ((4, 14), 0.625 * 100_000),
    ((15, 17), 0.8 * 100_000),
    ((18, 29), 0.86 * 100_000),
    ((30, 100), 0.95 * 100_000),
]:
    denmark_seroprevalences[age_range] = Prevalence(
        infections_per_100k=seroprevalence,
        date="1983",
        country="Denmark",
        active=Active.LATENT,
        source=DENMARK_SEROPREVALENCE_SOURCE,
    )

denmark_seroprevalence = Prevalence.weightedAverageByPopulation(
    (denmark_seroprevalences[0, 0], denmark_populations[0, 0]),
    (denmark_seroprevalences[1, 3], denmark_populations[1, 3]),
    (denmark_seroprevalences[4, 14], denmark_populations[4, 14]),
    (denmark_seroprevalences[15, 17], denmark_populations[15, 17]),
    (denmark_seroprevalences[18, 29], denmark_populations[18, 29]),
    (denmark_seroprevalences[30, 100], denmark_populations[30, 100]),
)


def estimate_prevalences():
    return [
        us_seroprevalence_2003_2010,
        denmark_seroprevalence,
        uk_seroprevalence_0_to_25,
    ]


def estimate_incidences():
    return []
