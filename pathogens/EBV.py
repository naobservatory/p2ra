from pathogen_properties import *

background = """EBV stands for Epstein-Barr virus. It is a type of 
herpes virus that infects humans and is known to cause infectious 
mononucleosis, also known as mono or glandular fever. 
EBV is a common virus that is transmitted through contact with infected
saliva, such as through kissing, sharing utensils, or close contact with
an infected person's respiratory droplets. EBV is a widespread virus that 
can persist in the body for life, although most people infected with EBV 
do not develop symptoms or have mild symptoms that resemble the flu."""

pathogen_chars = PathogenChars(
    na_type=NAType.DNA,
    enveloped=Enveloped.ENVELOPED,
    taxid=TaxID(10376),
)

denmark_population = Population(
    people=5_941_388,
    source="https://www.dst.dk/en/Statistik/emner/borgere/befolkning/befolkningstal",
    tag="Denmark-2023",
)

denmark_population_under_1 = Scalar(
    scalar=2 * 30_000,
    date="2023",
    source="https://www.dst.dk/en/Statistik/emner/borgere/befolkning/befolkningstal",
)

denmark_population_1_to_3 = Scalar(
    scalar=2 * 3 * 31_500,
    date="2023",
    source="https://www.dst.dk/en/Statistik/emner/borgere/befolkning/befolkningstal",
)

denmark_population_4_to_14 = Scalar(
    scalar=2 * 11 * 32_000,
    date="2023",
    source="https://www.dst.dk/en/Statistik/emner/borgere/befolkning/befolkningstal",
)

denmark_population_15_to_17 = Scalar(
    scalar=2 * 3 * 34_500,
    date="2023",
    source="https://www.dst.dk/en/Statistik/emner/borgere/befolkning/befolkningstal",
)

denmark_population_18_to_29 = Scalar(
    # Here, I actually added up each age rounded to the nearest 1000
    scalar=471_000 + 454_000,
    date="2023",
    source="https://www.dst.dk/en/Statistik/emner/borgere/befolkning/befolkningstal",
)

denmark_population_above_30 = Scalar(
    # I added all percentages for age groups 30-34 and up
    scalar=0.651 * denmark_population.people,
    date="2023",
    source="https://www.census.gov/popclock/world/da#:~:text=Population%20by%20Age%20and%20Sex",
)

denmark_seroprevalence = PrevalenceAbsolute(
    infections=(
        0.15 * denmark_population_under_1.scalar
        + 0.28 * denmark_population_1_to_3.scalar
        + 0.625 * denmark_population_4_to_14.scalar
        + 0.8 * denmark_population_15_to_17.scalar
        + 0.86 * denmark_population_18_to_29.scalar
        + 0.95 * denmark_population_above_30.scalar
    ),
    # This estimate applies a seroprevalence estimate from 1983 to Denmark's
    # current population.
    date="2023",
    active=Active.LATENT,
    tag="Denmark-2023",
    source="10.3109/inf.1983.15.issue-4.03",
)


def estimate_prevalences():
    return [denmark_seroprevalence.to_rate(denmark_population)]
