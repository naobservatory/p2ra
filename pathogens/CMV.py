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
overall_US_prevalence = Prevalence(
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


# Shedding duration data is from Table 3 of the paper linked below. I take all
# appplicable studies from the paper and then create a weighted average of
# their estimates.
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4494736/#:~:text=children%20and%20adults.-,Table%203,-Duration%20of%20CMV
adult_shedding_estimate = SheddingDuration(
    # The original article cited in the comment above estimates shedding
    # duration from figure 1C of the source. It says less than 90 days for all
    # patients. I estimate around 70 for the average patient. This is also
    # reasonably consistent with other inaccessible estimates from the paper
    # and other less reputable information online.
    days=70,
    number_of_participants=48,
    date=2004,
    source="https://pubmed.ncbi.nlm.nih.gov/2839977/#:~:text=The%20duration%20of%20CMV%20excretion%20varied%20from%203.0%20to%2028.4%20months%2C%20with%20a%20mean%20of%2013.0%20%2B/%2D%209.1%20months%20for%20urine",
)

child_shedding_estimate_1 = SheddingDuration(
    # Note that 30.4 is the # of days in an average month
    days=13 * 30.4,
    number_of_participants=79,
    date=1988,
    source="https://pubmed.ncbi.nlm.nih.gov/15529253/",
)

# I could not view the full paper, so I took information from the source paper
# linked in the comment above the shedding estimates
child_shedding_estimate_2 = SheddingDuration(
    # Note that 30.4 is the # of days in an average month
    days=6.5 * 30.4,
    number_of_participants=17,
    date=1970,
    source="https://pubmed.ncbi.nlm.nih.gov/4316198/",
)

child_shedding_average = SheddingDuration(
    days=(
        child_shedding_estimate_1.days
        * child_shedding_estimate_1.number_of_participants
        + child_shedding_estimate_2.days
        * child_shedding_estimate_2.number_of_participants
    )
    / (
        child_shedding_estimate_1.number_of_participants
        + child_shedding_estimate_2.number_of_participants
    ),
    number_of_participants=child_shedding_estimate_1.number_of_participants
    + child_shedding_estimate_2.number_of_participants,
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

# Seroprevalence of children ages 1-5
US_infant_seroprevalence = Prevalence(
    infections_per_100k=0.207 * 100_000,
    date="2011",
    active=Active.LATENT,
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8248283/#:~:text=the%20seroprevalence%20of%20cytomegalovirus%20immunoglobulin%20G%20(IgG)%20antibodies%20among%20US%20children%20aged%201%E2%80%935%20years%20was%2020.7%25",
)


US_child_prevalence = Prevalence(
    """To solve for the US_seroprevalence_of_children, we can simplify the 
equation US_child_proportion * US_child_prevalence + US_adult_proportion * 
new_US_adult_prevalence = overall_US_prevalence. 
This gives us that US_child_prevalence = (older_US_prevalence - 
US_adult_proportion * new_US_adult_prevalence) / US_child_proportion""",
    # This number evaluates to 0.29. The seroprevalence of children ages 1-5 was
    # 0.21, and it makes sense that this number is slightly lower than the overall
    # child prevalence, since most children are older than 5
    infections_per_100k=(
        overall_US_prevalence - US_adult_proportion * new_US_adult_prevalence
    )
    / US_child_proportion,
    date="2020",
    country="US",
    active=Active.LATENT,
    source="",
)

total_days_shedding_per_person = Scalar(
    scalar=US_child_proportion
    * US_child_prevalence
    * child_shedding_average.days
    + US_adult_proportion
    * new_US_adult_prevalence
    * adult_shedding_estimate.days,
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
    return [overall_US_prevalence, shedding_prevalence]
