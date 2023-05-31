from pathogen_properties import *

background = """Cytomegalovirus (CMV) is a common virus belonging to the 
herpesvirus family. It can infect people of all ages and is typically 
transmitted through close contact with infected body fluids, such as saliva, 
urine, blood, semen, and breast milk. CMV can also be transmitted during
pregnancy.

In healthy individuals with intact immune systems, CMV infection typically
lasts for a few weeks to a few months. However, the virus can remain dormant
in the body for life  after the acute phase, and may reactivate later under
certain circumstances such as when the immune system is weakened.
"""
pathogen_chars = PathogenChars(
    na_type=NAType.DNA,
    enveloped=Enveloped.ENVELOPED,
    taxid=TaxID(10359),
)

ger_adult_seroprevalence_estimate = Prevalence(
    # Estimate is for 18-79 year olds. In reality, this implies slightly lower
    # than 57% for the overall population including children
    infections_per_100k=0.567 * 100000,
    confidence_interval=(0.548 * 100000, 0.587 * 100000),
    number_of_participants=6552,
    country="Germany",
    date="1998",
    active=Active.LATENT,
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6059406/#:~:text=Overall%20CMV%20seroprevalence,Germany%20in%20women",
)

adult_seroprevalence_raleigh_durham_US = Prevalence(
    infections_per_100k=0.56 * 100_000,
    # Overall CMV seroprevalence was 56.7% (95%CI: 54.8–58.7%), with a higher seroprevalence in women (62.3%) than in men (51.0%). Seroprevalence increased with age: from 31.8% to 63.7% in men and from 44.1% to 77.6% in women when comparing the 18–29 with the 70–79 year age-group
    number_of_participants=694,
    date="2020",
    country="US",
    state="North Carolina",
    # Study is actually from the Raleigh-Durham-Chapel Hill Metropolitan Area
    county="Raleigh County",
    active=Active.LATENT,
    source="https://www.frontiersin.org/articles/10.3389/fcimb.2020.00334/full#:~:text=Methods%3A%20We,population%20was%2056%25",
)

nhanes_6_to_49_US_seroprevalence = Prevalence(
    infections_per_100k=0.504 * 100_000,
    country="US",
    active=Active.LATENT,
    number_of_participants=15_310,
    source="https://pubmed.ncbi.nlm.nih.gov/20426575/#:~:text=the%20overall%20age%2Dadjusted%20CMV%20seroprevalence%20was%2050.4%25",
    start_date="1999",
    end_date="2004",
)


def estimate_prevalences():
    return [
        nhanes_6_to_49_US_seroprevalence,
        ger_adult_seroprevalence_estimate,
        adult_seroprevalence_raleigh_durham_US,
    ]


def estimate_incidences():
    return []
