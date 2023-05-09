from pathogen_properties import *

background = """According to the CDC, “Rhinoviruses are the most frequent 
cause of the common cold. In the United States children have an average of two 
rhinovirus infections each year, and adults have an average of one.”
Furthermore, Rhinovirus declined somewhat less than other respiratory viruses
during the pandemic (https://jamanetwork.com/journals/jamanetworkopen
fullarticle/2801017).
Rhinovirus is not as seasonal as viruses like influenza or coronviruses.
(https://journals.plos.org/plosone/article/figure?id=10.1371/journal.pone.
0114871.g001)
(https://academic.oup.com/view-large/figure/89888454/195-6-773-fig003.jpeg)"""


pathogen_chars = PathogenChars(
    na_type=NAType.RNA,
    enveloped=Enveloped.NON_ENVELOPED,
    taxids=frozenset([TaxID(147711), TaxID(147712), TaxID(463676)]),
)

LA_county_population = Population(
    people=10_014_042,
    date="2020-04",
    tag="LA-2020",
    source="https://www.census.gov/quickfacts/fact/table/losangelescountycalifornia,CA/PST045221",
)


LA_county_under_18_population = Population(
    people=0.211 * LA_county_population.people,
    date="2020-04",
    tag="LA-2020",
    source="https://www.census.gov/quickfacts/fact/table/losangelescountycalifornia#:~:text=Persons%20under-,18,-years%2C%20percent",
)

LA_county_adult_population = Population(
    people=0.789 * LA_county_population.people,
    date="2020-04",
    tag="LA-2020",
    source="https://www.census.gov/quickfacts/fact/table/losangelescountycalifornia#:~:text=Persons%20under-,18,-years%2C%20percent",
)


"""This estimate was created by Simon with the following reasoning: During 
this study, "A single family respondent was contacted weekly by telephone to 
obtain information on the onset of acute respiratory or enteric illnesses in 
any of the household members. Across this time they tracked the number of 
respiratory illnesses, across age groups. These are overall respiratory 
illnesses (which we might equate with colds). In the same study, they report 
that colds were isolated in 34% of cases. (note that across 13140 cases of
respiratory illnesses, they only sampled 2227 of cases for isolation, which 
could lead to bias in this isolation rate). But for now, let's use this 
number. Checking a table of isolation rates across ages, the rate of 
Rhinoviruses among infections roughly holds true across ages."""
annual_rhinovirus_infections_under_18 = IncidenceAbsolute(
    # rhinovirus_share_of_cold = 0.34
    annual_infections=3.7 * 0.34 * LA_county_under_18_population.people,
    start_date="1976",
    end_date="1981",
    tag="LA-2020",
    source="doi.org/10.1017/S0950268800050779#?page=6",
)

annual_rhinovirus_infections_adults = IncidenceAbsolute(
    # rhinovirus_share_of_cold = 0.34
    annual_infections=2 * 0.34 * LA_county_adult_population.people,
    start_date="1976",
    end_date="1981",
    tag="LA-2020",
    source="doi.org/10.1017/S0950268800050779#?page=6",
)

rhinovirus_shedding_duration = SheddingDuration(
    days=12,
    confidence_interval=(10, 14),
    source="https://erj.ersjournals.com/content/44/1/169#:~:text=Virus%20shedding%20lasts%20on%20average%20for%2010%E2%80%9314%20days%20in%20immunocompetent%20subjects",
)

under_18_prevalence = annual_rhinovirus_infections_under_18.to_rate(
    LA_county_under_18_population
).to_prevalence(rhinovirus_shedding_duration)

adult_prevalence = annual_rhinovirus_infections_adults.to_rate(
    LA_county_adult_population
).to_prevalence(rhinovirus_shedding_duration)

total_prevalence = adult_prevalence + under_18_prevalence


"""This article analyzes information from the following studies, and their 
analysis seems to make sense from an initial screening. https://www.rcgp.org.
uk/clinical-and-research/our-programmes/research-and-surveillance-centre/
public-health-data#:~:text=December%202022%20%C2%A0(ZIP%20file%2C%203.7%20MB
"""
pandemic_decrease_factor = Scalar(
    scalar=0.1,
    date="2020",
    source="https://www.bmj.com/content/370/bmj.m3182#:~:text=still%20around%20nine%20times%20fewer%20cases%20than%20the%20five%20year%20average%20for%20this%20time%20of%20year.",
)


def estimate_prevalences():
    return [
        total_prevalence.__mul__(pandemic_decrease_factor),
    ]
