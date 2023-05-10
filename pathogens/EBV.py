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
    active=Active.LATENT,
    source="https://bmcpublichealth.biomedcentral.com/articles/10.1186/s12889-020-09049-x#:~:text=1982/2325%20individuals%20(85.3%25)%20were%20EBV%20seropositive",
)

# Children ages 8-19
nhanes_children_estimate = Prevalence(
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

US_child_proportion = Scalar(
    scalar=0.222,
    date="2022",
    source="https://www.census.gov/quickfacts/fact/table/US#",
)

US_adult_proportion = Scalar(
    scalar=1 - US_child_proportion.scalar,
    date="2022",
    source="https://www.census.gov/quickfacts/fact/table/US#",
)


# this estimate uses children 18-19 as a proxy for the adult population
us_seroprevalence_2003_2010 = nhanes_18_19_yo_estimate.__mul__(
    US_adult_proportion
).__add__(nhanes_children_estimate.__mul__(US_child_proportion))


def estimate_prevalences():
    return [us_seroprevalence_2003_2010]
