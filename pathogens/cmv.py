from pathogen_properties import *

background = """Cytomegalovirus (CMV) is a common virus belonging to the 
herpesvirus family. It can infect people of all ages and is typically 
transmitted through close contact with infected body fluids, such as saliva, 
urine, blood, semen, and breast milk. CMV can also be transmitted during
pregnancy.

Most people infected with CMV show no signs or symptoms. However, the virus 
remains dormant in latently infected cells and may reactivate later under 
certain circumstances, for instance when the immune system is weakened during 
organ transplantation.
"""
pathogen_chars = PathogenChars(
    na_type=NAType.DNA,
    enveloped=Enveloped.ENVELOPED,
    taxid=TaxID(10359),
)

de_adult_seroprevalence_estimate = Prevalence(
    # Estimate is for 18-79 year olds. In reality, this implies slightly lower
    # than 57% for the overall population including children
    infections_per_100k=0.567 * 100_000,
    confidence_interval=(0.548 * 100_000, 0.587 * 100_000),
    number_of_participants=6552,
    country="Germany",
    date="1998",
    active=Active.LATENT,
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6059406/#:~:text=Overall%20CMV%20seroprevalence,Germany%20in%20women",
)

adult_seroprevalence_raleigh_durham_US = Prevalence(
    infections_per_100k=0.56 * 100_000,
    # Overall CMV seroprevalence was 56.7% (95%CI: 54.8–58.7%), with a higher
    # seroprevalence in women (62.3%) than in men (51.0%). Seroprevalence
    # increased with age: from 31.8% to 63.7% in men and from 44.1% to 77.6%
    #  in women when comparing the 18–29 with the 70–79 year age-group
    number_of_participants=694,
    date="2020",
    country="United States",
    state="North Carolina",
    # Study is actually from the Raleigh-Durham-Chapel Hill Metropolitan Area
    county="Raleigh County",
    active=Active.LATENT,
    source="https://www.frontiersin.org/articles/10.3389/fcimb.2020.00334/full#:~:text=Methods%3A%20We,population%20was%2056%25",
)

nhanes_6_to_49_US_seroprevalence = Prevalence(
    infections_per_100k=0.504 * 100_000,
    country="United States",
    active=Active.LATENT,
    number_of_participants=15_310,
    source="https://pubmed.ncbi.nlm.nih.gov/20426575/#:~:text=the%20overall%20age%2Dadjusted%20CMV%20seroprevalence%20was%2050.4%25",
    start_date="1999",
    end_date="2004",
)

nhanes_6_to_49_US_urine_shedding = Prevalence(
    infections_per_100k=0.0383 * 100_000,
    # Among 6,828 CMV IgG-positive subjects tested, 537 had CMV DNA detected
    # in urine—a shedding prevalence of 9.70%. Among persons 6–49 years,
    # shedding prevalence was 3.83%. The prevalence of urinary shedding was
    # inversely associated with increasing age (26.60%, 6.50%, and 3.45% in
    # CMV IgG-positive subjects aged 6–11, 12–19, and 20–49 years, respectively
    country="United States",
    active=Active.ACTIVE,
    number_of_participants=6828,
    start_date="1999",
    end_date="2004",
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6097960/#:~:text=Among%206%2C828%20CMV%20IgG%2Dpositive%20subjects%20tested%2C%20537,11%2C%2012%E2%80%9319%2C%20and%2020%E2%80%9349%20years%2C%20respectively",
)


def estimate_prevalences():
    return [
        nhanes_6_to_49_US_seroprevalence,
        de_adult_seroprevalence_estimate,
        adult_seroprevalence_raleigh_durham_US,
    ]


def estimate_incidences():
    return []
