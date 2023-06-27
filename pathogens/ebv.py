import csv
import dataclasses
from collections import Counter

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
    selection=SelectionRound.ROUND_1,
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


CENSUS_QUERY = "https://data.census.gov/table?q=Annual+Estimates+of+the+Resident+Population+by+Single+Year"


race_cohorts: dict[str, dict[str, int]] = {
    "black_age_cohorts": {},
    "latino_age_cohorts": {},
    "white_age_cohorts": {},
}


digit_extractor = re.compile("\d+")

for race, cohort_dict in race_cohorts.items():
    with open(prevalence_data_filename(f"{race}.csv")) as inf:
        # Data is downloaded from this results page: "https://data.census.gov/table?q=Annual+Estimates+of+the+Resident+Population+by+Single+Year&g=010XX00US&tid=DECENNIALDHC2020.PCT12A".
        # Specifically, we used PCT12I (White alone, not Hispanic or Latino),
        # PCT12H (Hispanic or Latino), and PCT12B (Black or African American
        # Alone), saved as [white, latino, black]_age_cohorts.csv.
        reader = csv.reader(inf)

        # Skip first two non-age lines
        next(reader)
        next(reader)

        for row in reader:
            cohort_age = row[0]
            number = int(row[1].replace(",", ""))

            # Check if the row is a gender specification row
            if cohort_age.strip() in ["Male:", "Female:"]:
                continue

            # Set cohort_age to '0' for 'Under 1 year' label
            if cohort_age.strip() == "Under 1 year":
                cohort_age = "0"
            else:
                # Extract all integers and take the first one as age.
                age_digits = digit_extractor.findall(cohort_age)
                cohort_age = age_digits[0] if age_digits else None

            if cohort_age is not None:
                if int(cohort_age) <= 18:
                    if cohort_age in cohort_dict:
                        cohort_dict[cohort_age] += number
                    else:
                        cohort_dict[cohort_age] = number
                else:
                    if "adults" in cohort_dict:
                        cohort_dict["adults"] += number
                    else:
                        cohort_dict["adults"] = number


EBV_US_2003_2010_SEROPREVALENCE = "https://academic.oup.com/jid/article/208/8/1286/2192838#:~:text=Table%201.Demographic%20Factors%20Associated%20With%20Epstein%E2%80%93Barr%20Virus%20(EBV)%20Antibody%20(Ab)%20Prevalence%2C%20by%20Race/Ethnicity%E2%80%94National%20Health%20and%20Nutrition%20Examination%20Survey%20Cycles%202003%E2%80%932004%2C%202005%E2%80%932006%2C%202007%E2%80%932008%2C%20and%202009%E2%80%932010"


def us_seroprevalence_2020() -> Prevalence:
    with open(
        prevalence_data_filename("ebv_6_19_nhanes_2003_2010.csv")
    ) as inf:
        # Extracted from Table 1 of the above paper. This study aggregates NHANES
        # EBV seroprevalence measurements from 2003-2010 for 5 to 19 year olds.
        # We do not have seroprevalence data for 19+ year olds; given that EBV
        # stays latent after infection we treat the 18 to 19 year old cohort
        # prevalence as the prevalence of the overall adult US population.

        # We furthermore rescale the prevalence among different ethnic groups in
        # the US NHANES data by national shares of different ethnicities

        ethnicity_mapping = {
            "Mexican American": "latino_age_cohorts",
            "Black": "black_age_cohorts",
            "White": "white_age_cohorts",
        }
        estimate_weights: list[tuple[Prevalence, Population]] = []
        for row in csv.reader(inf):
            age_range, ethnicity, size, mean, ci_low, ci_high = [
                x.strip() for x in row
            ]
            if age_range in ["Age", "Total"]:
                continue

            prevalence = Prevalence(
                infections_per_100k=float(mean)
                * 1_000,  # percentage points to per
                # 100k,
                confidence_interval=(
                    float(ci_low) * 1_000,
                    float(ci_high) * 1_000,
                ),
                coverage_probability=0.95,
                number_of_participants=int(size),
                country="United States",
                start_date="2003",
                end_date="2010",
                active=Active.LATENT,
                source=EBV_US_2003_2010_SEROPREVALENCE,
            )

            if ethnicity in ethnicity_mapping:
                cohort_size = 0
                if (
                    age_range == "18_19"
                ):  # Matching 18-19 year old EBV+ rates with adult population
                    population = Population(
                        people=race_cohorts[ethnicity_mapping[ethnicity]][
                            "adults"
                        ],
                        source=CENSUS_QUERY,
                        country="United States",
                        date="2020",
                    )

                    continue

                population = Population(
                    people=sum(
                        race_cohorts[ethnicity_mapping[ethnicity]].values()
                    ),
                    source=CENSUS_QUERY,
                    country="United States",
                    date="2020",
                )

                estimate_weights.append((prevalence, population))

    return Prevalence.weightedAverageByPopulation(*estimate_weights)


def load_denmark_populations_by_age() -> dict[tuple[int, int], Population]:
    with open(prevalence_data_filename("DenmarkPopulationData.csv")) as file:
        # Source for CSV: https://www.census.gov/data-tools/demo/idb/#/pop?COUNTRY_YEAR=2023&COUNTRY_YR_ANIM=2023&FIPS_SINGLE=DA&menu=popViz&FIPS=DA&POP_YEARS=2023&popPages=BYAGE
        # We are applying the seroprevalence data from 1983 to current
        # population data, given that there is no more recent seroprevalence
        # data available.
        denmark_age_groups: Counter[tuple[int, int]] = Counter()
        # (min_age, max_age) -> people within that age range
        reader = csv.reader(file)
        next(reader)  # Skip the header

        for row in reader:
            # Skip the total population row
            if row[4] == "TOTAL":
                continue

            # Extract the age from the first column & turn "100+" into 100
            age = int(row[4].replace("+", ""))
            if age == 0:
                age_group = 0, 0
            elif age <= 3:
                age_group = 1, 3
            elif age <= 14:
                age_group = 4, 14
            elif age <= 17:
                age_group = 15, 17
            elif age <= 29:
                age_group = 18, 29
            elif age <= 100:
                age_group = 30, 100
            else:
                assert False
            # Remove spaces from the numbers and add them to the respective
            # age group
            denmark_age_groups[age_group] += int(row[5].replace(" ", ""))

        denmark_populations = {}
        for (min_age, max_age), people in denmark_age_groups.items():
            denmark_populations[min_age, max_age] = Population(
                people=people,
                date="2023",
                country="Denmark",
                source="https://www.census.gov/data-tools/demo/idb/#/pop?COUNTRY_YEAR=2023&COUNTRY_YR_ANIM=2023&FIPS_SINGLE=DA&menu=popViz&FIPS=DA&POP_YEARS=2023&popPages=BYAGE",
            )
    return denmark_populations


def denmark_seroprevalence_2023() -> Prevalence:
    DENMARK_SEROPREVALENCE_SOURCE = "10.3109/inf.1983.15.issue-4.03"

    denmark_seroprevalences: dict[tuple[int, int], Prevalence] = {}
    for (min_year, max_year), seroprevalence in [
        ((0, 0), 0.15 * 100_000),
        ((1, 3), 0.28 * 100_000),
        ((4, 14), 0.625 * 100_000),
        ((15, 17), 0.8 * 100_000),
        ((18, 29), 0.86 * 100_000),
        ((30, 100), 0.95 * 100_000),
    ]:
        denmark_seroprevalences[min_year, max_year] = Prevalence(
            infections_per_100k=seroprevalence,
            date="1983",
            country="Denmark",
            active=Active.LATENT,
            source=DENMARK_SEROPREVALENCE_SOURCE,
        )

    denmark_populations = load_denmark_populations_by_age()

    return Prevalence.weightedAverageByPopulation(
        *[
            (denmark_seroprevalences[age], denmark_populations[age])
            for age in denmark_populations.keys()
        ]
    )


def estimate_prevalences():
    denmark_2023 = denmark_seroprevalence_2023()
    # Seroprevalence should remain constant, so we can extrapolate from 1983
    # data, applied to the 2023 Denmark population backwards to 2021-2020.
    # This assumes that the population breakdown by age remained constant
    # between 2023 and 2021/2020.

    denmark_2020 = dataclasses.replace(
        denmark_2023, date_source=Variable(date="2020")
    )

    denmark_2021 = dataclasses.replace(
        denmark_2023, date_source=Variable(date="2021")
    )
    # Similar to Denmark, seroprevalence should remain constant, so we can
    # extrapolate from 2003-2010 data, applied to 2020 US population to 2021.
    us_2020 = us_seroprevalence_2020()
    us_2021 = dataclasses.replace(us_2020, date_source=Variable(date="2021"))

    return [
        us_2020,
        us_2021,
        denmark_2020,
        denmark_2021,
        denmark_2023,
        uk_seroprevalence_0_to_25,
    ]


def estimate_incidences():
    return []
