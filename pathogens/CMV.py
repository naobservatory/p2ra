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

# There is currently no vaccine for cytomegalovirus, so we can use
# seroprevalence to create a decent estimate of the overall proportion of
# people shedding at any given time.


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

# This study is from 2020, and finds a 56% seroprevalence for US adults.
adult_prevalence_raleigh_durham = Prevalence(
    infections_per_100k=0.56 * 100_000,
    number_of_participants=694,
    date="2020",
    country="US",
    state="North Carolina",
    # Study is actually from the Raleigh-Durham-Chapel Hill Metropolitan Area
    county="Raleigh County",
    active=Active.LATENT,
    source="https://www.frontiersin.org/articles/10.3389/fcimb.2020.00334/full#:~:text=Methods%3A%20We,population%20was%2056%25",
)

# This study is generally good, but also from around 20 years ago.
# It covers individuals ages 6-49.
# Since the Germany and North Carolina estimates are just for adults, and
# the NHANES US seroprevalence estimate covers a broader sample of the
#  population (and is pretty much in agreement with the other estimates), I
# think that the NHANES' 50% makes a lot of sense as an overall estimate.
nhanes_6_to_49_US_seroprevalence = Prevalence(
    infections_per_100k=0.504 * 100_000,
    country="US",
    active=Active.LATENT,
    number_of_participants=15_310,
    source="https://pubmed.ncbi.nlm.nih.gov/20426575/#:~:text=the%20overall%20age%2Dadjusted%20CMV%20seroprevalence%20was%2050.4%25",
    start_date="1999",
    end_date="2004",
)


# Shedding duration data is from Table 3 of the paper linked below. I take all
# appplicable studies from the paper and then create a weighted average of
# their estimates.
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4494736/#:~:text=children%20and%20adults.-,Table%203,-Duration%20of%20CMV
adult_shedding_estimate = SheddingDuration(
    # The original article cited in the comment above estimates shedding
    # duration from figure 1C of the source. It says less than 90 days for all
    # patients. I estimate around 60 for the average patient. This is also
    # reasonably consistent with other inaccessible estimates from the paper
    # linked above and other less reputable information online.
    days=60,
    number_of_participants=48,
    date="2004",
    source="https://academic.oup.com/jid/article/190/11/1908/836595#:~:text=Antibody%20maturation%20and%20presence%20of%20cytomegalovisrus%20(CMV)%20DNA%20in%20serum%2C%20after%20the%20onset%20of%20primary%20CMV%20infection%2C%20in%20immunocompetent%20patients%20(A%20and%20C)",
)

child_shedding_estimate_1 = SheddingDuration(
    # Note that 30.4 is the # of days in an average month
    days=13 * 30.4,
    number_of_participants=28,
    confidence_interval=(3.9 * 30.4, 22.1 * 30.4),  # standard deviation
    date="1988",
    source="https://pubmed.ncbi.nlm.nih.gov/2839977/#:~:text=The%20duration%20of%20CMV%20excretion%20varied%20from%203.0%20to%2028.4%20months%2C%20with%20a%20mean%20of%2013.0",
)

# I could not view the full paper, so I took information from the source paper
# linked in the comment above the shedding estimates
child_shedding_estimate_2 = SheddingDuration(
    # Note that 30.4 is the # of days in an average month
    # In parentheses is the number of months the average participant shed
    days=(13 / 17 * 10.5 + 2 / 17 * 12 + 1 / 17 * 18 + 1 / 17 * 19) * 30.4,
    number_of_participants=17,
    date="1970",
    # See the 3rd paragraph of page 412 for relevant shedding information
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
# This estimate is not part of our final estimate, but is used to fact check that the child prevalence estimate below is reasonable
US_1_to_5_seroprevalence = Prevalence(
    infections_per_100k=0.207 * 100_000,
    date="2011",
    number_of_participants=698,
    active=Active.LATENT,
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8248283/#:~:text=the%20seroprevalence%20of%20cytomegalovirus%20immunoglobulin%20G%20(IgG)%20antibodies%20among%20US%20children%20aged%201%E2%80%935%20years%20was%2020.7%25",
)


# To solve for the US_child_prevalence, we can first look at how to
# calculate overall US prevalence:
# nhanes_6_to_49_US_seroprevalence (taken to be overall prevalence) =
# US_child_proportion * US_child_prevalence + US_adult_proportion
# * adult_prevalence_raleigh_durham
# Based on this, we can solve for US_child_prevalence the following way:
#  US_child_prevalence = (nhanes_6_to_49_US_seroprevalence -  US_adult_proportion * adult_prevalence_raleigh_durham)) /
# US_child_proportion.
# This number evaluates to 0.29. The seroprevalence of children ages 1-5
# was 0.21. It makes sense that this number is slightly lower than the
# overall child prevalence, since most children are older than 5.
# Note that we have slightly less confidence in
# adult_prevalence_releigh_durham, because it is not meant to fit the
# entire population like the NHANES estimate is
US_child_prevalence = (
    (
        nhanes_6_to_49_US_seroprevalence
        - (adult_prevalence_raleigh_durham * US_adult_proportion)
    )
    / (US_child_proportion)
).target(date="2022")


total_days_shedding_per_person = Scalar(
    scalar=US_child_proportion.scalar
    * US_child_prevalence.infections_per_100k
    * child_shedding_average.days
    + US_adult_proportion.scalar
    * adult_prevalence_raleigh_durham.infections_per_100k
    * adult_shedding_estimate.days,
    source="",
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
    return [nhanes_6_to_49_US_seroprevalence, shedding_prevalence]
