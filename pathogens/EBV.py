import csv
from collections import Counter

from pathogen_properties import *

background = """EBV stands for Epstein-Barr virus. It is a type of 
herpes virus that infects humans and is known to cause infectious 
mononucleosis, also known as mono or glandular fever. 
EBV is a common virus that is transmitted through contact with infected
saliva, such as through kissing, sharing utensils, or close contact with
an infected person's respiratory droplets. EBV is a widespread virus that 
can persist in the body for life, although most people infected with EBV 
do not develop symptoms or have mild symptoms that resemble the flu."""

# TODO: Simon will look into incidence of mononucleosis in the US. There is
# also the question of mono being a proxy for EBV infection, or activation.

pathogen_chars = PathogenChars(
    na_type=NAType.DNA,
    enveloped=Enveloped.ENVELOPED,
    taxid=TaxID(10376),
)

# "We conducted a sero-epidemiological survey using serum samples from 2325
# individuals between 0 and 25 years old to assess prevalence of detectable
# anti-EBV antibodies."
UK_seroprevalence_0_to_25 = Prevalence(
    # This study is not used in the estimate,
    # but is used to check that the later estimate is reasonable
    infections_per_100k=0.853 * 100_000,
    country="UK",
    start_date="2002",
    end_date="2013",
    number_of_participants=2325,
    active=Active.LATENT,
    source="https://bmcpublichealth.biomedcentral.com/articles/10.1186/s12889-020-09049-x#:~:text=1982/2325%20individuals%20(85.3%25)%20were%20EBV%20seropositive",
)

# Children ages 6-19
nhanes_age_6_19_estimate = Prevalence(
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

nhanes_18_19_yo_estimate = Prevalence(
    # This number comes from 2009-2010 and pulls only the most recent data
    # from the Balfour NHANES estimate above
    infections_per_100k=0.89 * 100_000,
    start_date="2003",
    end_date="2010",
    country="United States",
    active=Active.LATENT,
    source="https://academic.oup.com/jid/article/208/8/1286/2192838#:~:text=years%2C%2069%25%3B%20and-,18%E2%80%9319%20years%2C%2089%25,-.%20Within%20each%20race",
)

under_18_population_US = Population(
    people=0.222 * 333_287_557,
    date="2022",
    tag="US-2022",
    source="https://www.census.gov/quickfacts/fact/table/US#",
)

over_18_population_US = Population(
    people=(333_287_557 - under_18_population_US.people),
    date="2022",
    tag="US-2022",
    source="https://www.census.gov/quickfacts/fact/table/US#",
)

# this estimate uses children 18-19 as a proxy for the adult population
us_seroprevalence_2003_2010 = Prevalence.weightedAverageByPopulation(
    (nhanes_age_6_19_estimate, under_18_population_US),
    (nhanes_18_19_yo_estimate, over_18_population_US),
)

denmark_population = Population(
    people=5_941_388,
    source="https://www.dst.dk/en/Statistik/emner/borgere/befolkning/befolkningstal",
    tag="Denmark-2023",
)

denmark_age_groups: Counter[tuple[int, int]] = Counter()

# Source for CSV: https://www.census.gov/data-tools/demo/idb/#/pop?COUNTRY_YEAR=2023&COUNTRY_YR_ANIM=2023&FIPS_SINGLE=DA&menu=popViz&FIPS=DA&POP_YEARS=2023&popPages=BYAGE
with open(prevalence_data_filename("DenmarkPopulationData.csv")) as file:
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
        population = float(row[5].replace(" ", ""))

        if age == 0:
            denmark_age_groups[(0, 0)] += int(population)
        elif 1 <= age <= 3:
            denmark_age_groups[(1, 3)] += int(population)
        elif 4 <= age <= 14:
            denmark_age_groups[(4, 14)] += int(population)
        elif 15 <= age <= 17:
            denmark_age_groups[(15, 17)] += int(population)
        elif 18 <= age <= 29:
            denmark_age_groups[(18, 29)] += int(population)
        elif 30 <= age <= 100:
            denmark_age_groups[(30, 100)] += int(population)
        else:
            assert False

    denmark_populations = {}
    for (min_age, max_age), people in denmark_age_groups.items():
        denmark_populations[min_age, max_age] = Population(
            people=people,
            date="2023",
            country="Denmark",
            tag="%s-%syo" % (min_age, max_age),
            source="https://www.census.gov/data-tools/demo/idb/#/pop?COUNTRY_YEAR=2023&COUNTRY_YR_ANIM=2023&FIPS_SINGLE=DA&menu=popViz&FIPS=DA&POP_YEARS=2023&popPages=BYAGE",
        )
    print(denmark_populations[0, 0])


# These estimates apply seroprevalence estimates from 1983 to Denmark's
# current population.
denmark_seroprevalences = {
    (0, 0): Prevalence(
        infections_per_100k=0.15 * 100_000,
        date="2023",
        country="Denmark",
        active=Active.LATENT,
        # See page 336 of the study (page 3 of the pdf)
        source="10.3109/inf.1983.15.issue-4.03",
    ),
    (1, 3): Prevalence(
        infections_per_100k=0.28 * 100_000,
        date="2023",
        country="Denmark",
        active=Active.LATENT,
        # See page 336 of the study (page 3 of the pdf)
        source="10.3109/inf.1983.15.issue-4.03",
    ),
    (4, 14): Prevalence(
        infections_per_100k=0.625 * 100_000,
        date="2023",
        country="Denmark",
        active=Active.LATENT,
        # See page 336 of the study (page 3 of the pdf)
        source="10.3109/inf.1983.15.issue-4.03",
    ),
    (15, 17): Prevalence(
        infections_per_100k=0.8 * 100_000,
        date="2023",
        country="Denmark",
        active=Active.LATENT,
        # See page 336 of the study (page 3 of the pdf)
        source="10.3109/inf.1983.15.issue-4.03",
    ),
    (18, 29): Prevalence(
        infections_per_100k=0.86 * 100_000,
        date="2023",
        country="Denmark",
        active=Active.LATENT,
        # See page 336 of the study (page 3 of the pdf)
        source="10.3109/inf.1983.15.issue-4.03",
    ),
    (30, 100): Prevalence(
        infections_per_100k=0.95 * 100_000,
        date="2023",
        country="Denmark",
        active=Active.LATENT,
        # See page 336 of the study (page 3 of the pdf)
        source="10.3109/inf.1983.15.issue-4.03",
    ),
}

denmark_overall_seroprevalence = Prevalence.weightedAverageByPopulation(
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
        denmark_overall_seroprevalence,
        UK_seroprevalence_0_to_25,
    ]
