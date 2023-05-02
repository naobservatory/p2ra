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

UK_study = Prevalence(
    infections_per_100k=85300,
    country="UK",
    date="2002-2013",
    active=Active.LATENT(),
)

# this study cites a  textbook published in 2007. It is not used
# in the estimate, but rather to confirm that the estimate is reasonable
global_estimate = Prevalence(
    infections_per_100k=90000,
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4235690/#:~:text=It%20is%20estimated%20that%20over%2090%25%20of%20the%20world’s%20population%20is%20infected%20with%20the%20virus",
)

# this study looks only at children. It is also not used in the estimate,
# but helps confirm that the estimate is reasonable
US_children_study = Prevalence(
    infections_per_100k=66500,
    country="United States",
    date="2003-2010",
    source="https://pubmed.ncbi.nlm.nih.gov/23717674/#:~:text=Overall%20EBV%20seroprevalence%20was%2066.5%25%20(95%25%20CI%2064.3%25%2D68.7%25.)",
)


def estimate_prevalences():
    return [UK_study]
