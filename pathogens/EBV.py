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
    infections_per_100k=0.853 * 100_000,
    country="UK",
    start_date="2002",
    end_date="2013",
    active=Active.LATENT,
    source="https://bmcpublichealth.biomedcentral.com/articles/10.1186/s12889-020-09049-x#:~:text=1982/2325%20individuals%20(85.3%25)%20were%20EBV%20seropositive",
)

# this study cites a  textbook published in 2007. It is not used
# in the estimate, but rather to confirm that the estimate is reasonable
global_estimate = Prevalence(
    infections_per_100k=0.9 * 100_000,
    active=Active.LATENT,
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4235690/#:~:text=It%20is%20estimated%20that%20over%2090%25%20of%20the%20world’s%20population%20is%20infected%20with%20the%20virus",
)


nhanes_children_estimate = Prevalence(
    # This study is also not used in the estimate,
    # but helps confirm that the estimate is reasonable
    infections_per_100k=0.665 * 100_000,
    confidence_interval=(0.643 - 0.687),
    country="United States",
    start_date="2003",
    end_date="2010",
    number_of_participants=8417,
    active=Active.LATENT,
    source="https://pubmed.ncbi.nlm.nih.gov/23717674/#:~:text=Overall%20EBV%20seroprevalence%20was%2066.5%25%20(95%25%20CI%2064.3%25%2D68.7%25.)",
)


def estimate_prevalences():
    return [nhanes_children_estimate]
