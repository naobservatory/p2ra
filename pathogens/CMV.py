from pathogen_properties import *

background = """Cytomegalovirus (CMV) is a common virus belonging to the 
herpesvirus family. It can infect people of all ages and is typically 
transmitted through close contact with infected body fluids, such as saliva, 
urine, blood, semen, and breast milk. CMV can also be transmitted from a 
pregnant person to their fetus during pregnancy.

In healthy individuals with intact immune systems, CMV infection typically
lasts for a few weeks to a few months. The acute phase of the infection, 
during which the virus replicates and causes symptoms, usually lasts for 
several weeks. However, the virus can remain dormant in the body for life 
after the acute phase, and may reactivate later under certain circumstances
such as when the immune system is weakened. According to Mount Sinai Hospital, 
cytomegalovirus lasts for around 5 weeks. According to the CDC, about 1 in 200 
babies are born with cytomegalovirus.  
"""
pathogen_chars = PathogenChars(
    na_type=NAType.DNA,
    enveloped=Enveloped.ENVELOPED,
    taxid=TaxID(10359),
)

# There is currently no vaccine for cytomegalovirus, so we can use seroprevalence to create a decent estimate of the overall proportion of people shedding at any given time.


GER_adult_seroprevalence_estimate = Prevalence(
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

# This study is generally good, but also from around 20 years ago
older_US_prevalence = Prevalence(
    infections_per_100k=0.504 * 100_000,
    country="US",
    active=Active.LATENT,
    source="https://pubmed.ncbi.nlm.nih.gov/20426575/#:~:text=the%20overall%20age%2Dadjusted%20CMV%20seroprevalence%20was%2050.4%25",
    start_date="1999",
    end_date="2004",
)

# This study is from 2020, and finds a 56% seroprevalence for US adults.
new_US_adult_prevalence = Prevalence(
    infections_per_100k=0.56 * 100_000,
    date="2020",
    country="US",
    active=Active.LATENT,
    source="https://www.frontiersin.org/articles/10.3389/fcimb.2020.00334/full#:~:text=Methods%3A%20We,population%20was%2056%25",
)

# Since the 2 adult seroprevalence estimates are around 56-57%, these
# estimates are slightly too high because they don't account for children,
# and the 1 overall seroprevalence estimate is 50%, I think that 50% makes
# a lot of sense as an overall estimate


# data is from Table 3 of the source paper. This is not an exact number, but an estimate based on the multiple estimates cited in the paper (and the paper's own analysis)
adult_shedding = SheddingDuration(
    days=90,
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4494736/#:~:text=children%20and%20adults.-,Table%203,-Duration%20of%20CMV",
)


# data is from Table 3 of the source paper. This is not an exact number, but an estimate based on the multiple estimates cited in the paper (and the paper's own analysis)
child_shedding = SheddingDuration(
    days=335,
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4494736/#:~:text=children%20and%20adults.-,Table%203,-Duration%20of%20CMV",
)

"""The CDC says that in the United States, nearly one in three children is
already infected with CMV by age five.
Over half of adults have been infected with CMV by age 40. I estimate that 40% 
of the population gets CMV as a child, so according to our seroprevalence 
estimate, 57% - 40% = 17% gets CMV as an adult"""
total_days_shedding_per_person = Scalar(
    scalar=0.40 * child_shedding.days + 0.17 * adult_shedding.days,
    source="https://www.cdc.gov/cmv/overview.html#:~:text=nearly%20one%20in%20three%20children%20is%20already%20infected%20with%20CMV%20by%20age%20five.%20Over%20half%20of%20adults%20have%20been%20infected%20with%20CMV%20by%20age%2040.",
)

# average lifespan in the US in days - average lifespan in years * days/year
average_lifespan_US = Scalar(
    scalar=76.1 * 365,
    source="https://www.cdc.gov/nchs/pressroom/nchs_press_releases/2022/20220831.htm",
)

shedding_prevalence = Prevalence(
    infections_per_100k=total_days_shedding_per_person.scalar
    / average_lifespan_US.scalar
    * 100000,
    active=Active.ACTIVE,
)


def estimate_prevalences():
    return [GER_adult_seroprevalence_estimate, shedding_prevalence]
