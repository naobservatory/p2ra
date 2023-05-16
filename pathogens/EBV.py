import csv

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

# Children ages 8-19
nhanes_age_8_19_estimate = Prevalence(
    infections_per_100k=0.665 * 100_000,
    confidence_interval=(0.643 * 100_000, 0.687 * 100_000),
    country="United States",
    start_date="2003",
    end_date="2010",
    number_of_participants=8417,
    active=Active.LATENT,
    source="https://pubmed.ncbi.nlm.nih.gov/23717674/#:~:text=Overall%20EBV%20seroprevalence%20was%2066.5%25%20(95%25%20CI%2064.3%25%2D68.7%25.)",
)

nhanes_18_19_yo_estimate = Prevalence(
    infections_per_100k=0.89 * 100_000,
    start_date="2003",
    end_date="2010",
    country="United States",
    active=Active.LATENT,
    source="https://academic.oup.com/jid/article/208/8/1286/2192838#:~:text=years%2C%2069%25%3B%20and-,18%E2%80%9319%20years%2C%2089%25,-.%20Within%20each%20race",
)

under_18_proportion_US = Scalar(
    scalar=0.222,
    date="2022",
    source="https://www.census.gov/quickfacts/fact/table/US#",
)

over_18_proportion_US = Scalar(
    scalar=1 - under_18_proportion_US.scalar,
    date="2022",
    source="https://www.census.gov/quickfacts/fact/table/US#",
)


# this estimate uses children 18-19 as a proxy for the adult population
us_seroprevalence_2003_2010 = nhanes_18_19_yo_estimate * (
    over_18_proportion_US
) + nhanes_age_8_19_estimate * (under_18_proportion_US)

denmark_population = Population(
    people=5_941_388,
    source="https://www.dst.dk/en/Statistik/emner/borgere/befolkning/befolkningstal",
    tag="Denmark-2023",
)

age_groups = {
    "denmark_population_under_1": Population(
        people=0,
        date="2023",
        country="Denmark",
        tag="0-1yo",
    ),
    "denmark_population_1_to_3": Population(
        people=0,
        date="2023",
        country="Denmark",
        tag="1-3yo",
    ),
    "denmark_population_4_to_14": Population(
        people=0,
        date="2023",
        country="Denmark",
        tag="4-14yo",
    ),
    "denmark_population_15_to_17": Population(
        people=0,
        date="2023",
        country="Denmark",
        tag="15-17yo",
    ),
    "denmark_population_18_to_29": Population(
        people=0,
        date="2023",
        country="Denmark",
        tag="18-29yo",
    ),
    "denmark_population_over_30": Population(
        people=0,
        date="2023",
        country="Denmark",
        tag="30+yo",
    ),
}

# Open the CSV file
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
        population = int(row[5].replace(" ", ""))

        if age == 0:
            age_groups["denmark_population_under_1"].people += population
        elif 1 <= age <= 3:
            age_groups["denmark_population_1_to_3"].people += population
        elif 4 <= age <= 14:
            age_groups["denmark_population_4_to_14"].people += population
        elif 15 <= age <= 17:
            age_groups["denmark_population_15_to_17"].people += population
        elif 18 <= age <= 29:
            age_groups["denmark_population_18_to_29"].people += population
        elif 30 <= age <= 100:
            age_groups["denmark_population_over_30"].people += population


for group, population in age_groups.items():
    print(f"{group}: {population}")

denmark_seroprevalence = PrevalenceAbsolute(
    infections=(
        0.15 * age_groups["denmark_population_under_1"].people
        + 0.28 * age_groups["denmark_population_1_to_3"].people
        + 0.625 * age_groups["denmark_population_4_to_14"].people
        + 0.8 * age_groups["denmark_population_15_to_17"].people
        + 0.86 * age_groups["denmark_population_18_to_29"].people
        + 0.95 * age_groups["denmark_population_above_30"].people
    ),
    # This estimate applies a seroprevalence estimate from 1983 to Denmark's
    # current population.
    date="2023",
    country="Denmark",
    active=Active.LATENT,
    tag="Denmark-2023",
    # See page 336 of the study (page 3 of the pdf)
    source="10.3109/inf.1983.15.issue-4.03",
)


def estimate_prevalences():
    return [
        us_seroprevalence_2003_2010,
        denmark_seroprevalence.to_rate(denmark_population),
    ]
