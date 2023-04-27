from pathogen_properties import *

background = """Cytomegalovirus (CMV) is a common virus belonging to the herpesvirus family. It can infect people of all ages and is typically transmitted through close contact with infected body fluids, such as saliva, urine, blood, semen, and breast milk. CMV can also be transmitted from a pregnant person to their fetus during pregnancy.

In healthy individuals with intact immune systems, CMV infection typically lasts for a few weeks to a few months. The acute phase of the infection, during which the virus replicates and causes symptoms, usually lasts for several weeks. However, the virus can remain dormant in the body for life after the acute phase, and may reactivate later under certain circumstances, such as when the immune system is weakened.

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


adult_shedding = Scalar(
    scalar=90,
    source="https://pubmed.ncbi.nlm.nih.gov/15529253",
)


def estimate_prevalences():
    return [Germany_prevalence, Germany_prevalence]
