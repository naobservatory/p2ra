from pathogen_properties import *

background = """Cytomegalovirus (CMV) is a common virus belonging to the 
herpesvirus family. It can infect people of all ages and is typically 
transmitted through close contact with infected body fluids, such as saliva, 
urine, blood, semen, and breast milk. CMV can also be transmitted from a pregnant person to their fetus during pregnancy.

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

# Studies patients 18-79 years old
GER_seroprevalence_estimate = Prevalence(
    infections_per_100k=0.567 * 100000,
    confidence_interval=(0.548 * 100000, 0.587 * 100000),
    number_of_participants=6552,
    country="Germany",
    date="1998",
    active=Active.LATENT,
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6059406/#:~:text=Overall%20CMV%20seroprevalence,Germany%20in%20women",
)

US_prevalence = Prevalence(
    # According to the CDC,  "over 50% of people in the US have been infected
    # with CMV by age 40." At age 5 1/3 of people are infected, and by age 40
    # 1/2 are infected, which hints toward the prevalence being asymptotic at
    # around 50%. This implies that the true proportion is somewhere between
    # around 1/3 and 3/5 of people. This number is not actually used, but can
    # help confirm that the Germany estimate is reasonable
    infections_per_100k=0.5 * 100_000,
    country="US",
    active=Active.LATENT,
    # The CDC's source on this information is unclear.
    source="https://www.cdc.gov/cmv/overview.html#:~:text=In%20the%20United%20States%2C%20nearly%20one%20in%20three%20children%20is%20already%20infected%20with%20CMV%20by%20age%20five.%20Over%20half%20of%20adults%20have%20been%20infected%20with%20CMV%20by%20age%2040.%20Once",
    date="2020",
)

# Since the US-specific number is very rough, I think it makes sense to take the Germany estimate of 57000 per 100k.


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
    source="https://www.cdc.gov/cmv/overview.html",
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
    return [GER_seroprevalence_estimate, shedding_prevalence]
