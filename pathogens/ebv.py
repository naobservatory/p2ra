import csv
from collections import Counter
from typing import Dict, List

import numpy as np

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


EBV_US_2003_2010_SEROPREVALENCE = "https://academic.oup.com/jid/article/208/8/1286/2192838#:~:text=Table%201.Demographic%20Factors%20Associated%20With%20Epstein%E2%80%93Barr%20Virus%20(EBV)%20Antibody%20(Ab)%20Prevalence%2C%20by%20Race/Ethnicity%E2%80%94National%20Health%20and%20Nutrition%20Examination%20Survey%20Cycles%202003%E2%80%932004%2C%202005%E2%80%932006%2C%202007%E2%80%932008%2C%20and%202009%E2%80%932010"

with open(prevalence_data_filename("ebv_6_19_nhanes_2003_2010.csv")) as inf:
    # Extracted from Table 1 of the above paper. This study aggregates NHANES
    # EBV seroprevalence measurements from 2003-2010 for 5 to 19 year olds.
    # We do not have seroprevalence data for 19+ year olds; given that EBV
    # stays latent after infection we treat the 18 to 19 year old cohort
    # prevalence as the prevalence of the overall adult US population.

    cols = None

    age_6_17: Dict[str, List] = {
        "mean": [],
        "ci_lower": [],
        "ci_upper": [],
        "weights": [],
        "cohort_sizes": [],
    }

    age_18_19: Dict[str, List] = {
        "mean": [],
        "ci_lower": [],
        "ci_upper": [],
        "weights": [],
        "cohort_sizes": [],
    }

    total_participants = 0

    for row in csv.reader(inf):
        if cols is None:
            cols = row
            continue
        age = row[0].strip()
        cohort_size = int(row[2])
        mean = float(row[3]) / 100
        ci_lower = float(row[4]) / 100
        ci_upper = float(row[5]) / 100

        if age == "Total":
            total_participants += cohort_size
            continue

        elif age != "18_19":
            age_6_17["mean"].append(mean)
            age_6_17["ci_lower"].append(ci_lower)
            age_6_17["ci_upper"].append(ci_upper)
            age_6_17["weights"].append(cohort_size / total_participants)
            age_6_17["cohort_sizes"].append(cohort_size)
            continue

        elif age == "18_19":
            age_18_19["mean"].append(mean)
            age_18_19["ci_lower"].append(ci_lower)
            age_18_19["ci_upper"].append(ci_upper)
            age_18_19["weights"].append(cohort_size / total_participants)
            age_18_19["cohort_sizes"].append(cohort_size)

    mean_18_19 = np.average(age_18_19["mean"], weights=age_18_19["weights"])
    ci_upper_18_19 = np.average(
        age_18_19["ci_upper"], weights=age_18_19["weights"]
    )
    ci_lower_18_19 = np.average(
        age_18_19["ci_lower"], weights=age_18_19["weights"]
    )

    mean_6_17 = np.average(age_6_17["mean"], weights=age_6_17["weights"])
    ci_upper_6_17 = np.average(
        age_6_17["ci_upper"], weights=age_6_17["weights"]
    )
    ci_lower_6_17 = np.average(
        age_6_17["ci_lower"], weights=age_6_17["weights"]
    )

    nhanes_6_17_yo_seroprevalence_2003_2010 = Prevalence(
        infections_per_100k=mean_6_17 * 100_000,
        confidence_interval=(ci_lower_6_17 * 100_000, ci_upper_6_17 * 100_000),
        coverage_probability=0.95,
        number_of_participants=sum(age_6_17["cohort_sizes"]),
        country="United States",
        start_date="2003",
        end_date="2010",
        active=Active.LATENT,
        source=EBV_US_2003_2010_SEROPREVALENCE,
    )

    nhanes_18_19_yo_seroprevalence_2003_2010 = Prevalence(
        infections_per_100k=mean_18_19 * 100_000,
        confidence_interval=(
            ci_lower_18_19 * 100_000,
            ci_upper_18_19 * 100_000,
        ),
        coverage_probability=0.95,
        number_of_participants=sum(age_18_19["cohort_sizes"]),
        country="United States",
        start_date="2003",
        end_date="2010",
        active=Active.LATENT,
        source=EBV_US_2003_2010_SEROPREVALENCE,
    )

us_seroprevalence_2003_2010 = Prevalence.weightedAverageByPopulation(
    (nhanes_6_17_yo_seroprevalence_2003_2010, us_population_u18),
    (nhanes_18_19_yo_seroprevalence_2003_2010, us_population_18plus),
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
