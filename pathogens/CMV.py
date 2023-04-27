from pathogen_properties import *

background = """Cytomegalovirus (CMV) is a common virus belonging to the herpesvirus family. 
It can infect people of all ages and is typically transmitted through close contact with infected body fluids,
 such as saliva, urine, blood, semen, and breast milk. CMV can also be transmitted from a pregnant person to their fetus during pregnancy.

In healthy individuals with intact immune systems, CMV infection typically lasts for a few weeks to a few months. 
The acute phase of the infection, during which the virus replicates and causes symptoms, usually lasts for several weeks. 
However, the virus can remain dormant in the body for life after the acute phase, and may reactivate later under certain circumstances,
such as when the immune system is weakened.
According to Mount Sinai Hospital, cytomegalovirus lasts for around 5 weeks.
According to the CDC, about half of babies are born with cytomegalovirus.  
"""

# There is currently no vaccine for cytomegalovirus, so we can use seroprevalence to create a decent estimate of the overall proportion of people shedding at any given time.

pathogen_chars = PathogenChars(
    na_type=NAType.DNA,
    enveloped=Enveloped.ENVELOPED,
    taxid=TaxID(10359),
)

# Studies patients 18-79 years old
Germany_prevalence = Prevalence(
    infections_per_100k=57000,
    country="Germany",
    date="2002-2013",
    active=Active.LATENT,
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6059406/",
)

# Specifically, "over 50% of people in the US have been infected with CMV by age 40." 40 is around the average age in the US, so this implies a prevalence of around 50%
US_prevalence = Prevalence(
    infections_per_100k=50000,
    country="US",
    source="https://www.cdc.gov/cmv/overview.html",
    date="Current",
)

US_Europe_Prevalence = Prevalence(
    infections_per_100k=90000,
    source="https://pubmed.ncbi.nlm.nih.gov/20564615/",
    date="2010",
)

# Since the US-specific number is non-exact and lower, and the US / Europe estimate is much higher and includes a very different set of countries, I think it makes sense to take the Germany estimate of 57000 per 100k


# how many days infected adults shed for on average
# data is from Table 3 of the source paper. This is not an exact number, but an estimate based on the multiple estimates cited in the paper (and the paper's own analysis)
adult_shedding = Scalar(
    scalar=90,
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4494736/",
)

# how many days infected children shed for on average
# data is from Table 3 of the source paper. This is not an exact number, but an estimate based on the multiple estimates cited in the paper (and the paper's own analysis)
children_shedding = Scalar(
    scalar=335,
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4494736/",
)

"""The CDC says that in the United States, nearly one in three children is already infected with CMV by age five.
Over half of adults have been infected with CMV by age 40
I estimate that 40% of the population gets CMV as a child, so according to our seroprevalence estimate, 17% gets CMV as an adult"""
total_days_shedding_per_person = Scalar(
    scalar=0.40 * children_shedding.scalar + 0.17 * adult_shedding.scalar,
    source="https://www.cdc.gov/cmv/overview.html",
)

# average lifespan in the US in days
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
    return [Germany_prevalence, shedding_prevalence]
