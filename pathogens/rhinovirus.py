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

# TODOs for SIMON:
#  - Check if national estimates can be applied to regions like Ohio during some specific time period.
#  - Check if not accounting for seasonality is reasonable

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

# Ideally, we would find the under 19 population instead to match our under 19
# prevalence estimate
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
annual_rhinovirus_infections_under_19 = IncidenceAbsolute(
    # The weighted average annual number of respiratory illnesses among
    # 0-19-year-olds is in parentheses, aggregating the mean of 0-4 and 5-19
    # year olds from table 1 of the source
    # rhinovirus_share_of_cold = 0.34
    annual_infections=((539 * 4.9 + 1541 * 2.8) / (539 + 1541))
    * 0.34
    * LA_county_under_18_population.people,
    start_date="1976",
    end_date="1981",
    tag="LA-2020",
    source="doi.org/10.1017/S0950268800050779#?page=6",
)

annual_rhinovirus_infections_adults = IncidenceAbsolute(
    # The weighted average annual number of respiratory illnesses among
    # 20+ year-olds is in parentheses, aggregating the mean of 20-40 and 40+
    # year olds from table 1 of the source
    # rhinovirus_share_of_cold = 0.34
    annual_infections=((1523 * 2.2 + 1757 * 1.6) / (1523 + 1757))
    * 0.34
    * LA_county_adult_population.people,
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

under_18_prevalence = annual_rhinovirus_infections_under_19.to_rate(
    LA_county_under_18_population
).to_prevalence(rhinovirus_shedding_duration)

adult_prevalence = annual_rhinovirus_infections_adults.to_rate(
    LA_county_adult_population
).to_prevalence(rhinovirus_shedding_duration)

total_prevalence = adult_prevalence + under_18_prevalence


# Link for all RCGP reports (specific reference is given under source):
#  https://www.rcgp.org.uk/clinical-and-research/our-programmes/research-and-surveillance-centre/public-health-data
pandemic_decrease_factor = Scalar(
    scalar=0.33,
    date="2020",
    # See page 8 of the December 2020 report, around weeks 44 to 53
    # The 5-year avg at this time is in blue, and the national avg in red below
    # From this graph, the UK's decrease factor was around 0.2.
    # On page 5, we see that Covid incidence in the relevant weeks
    # is around 14-28 per 100k per day.
    # During fall 2020, LA County had an incidence of around 10-150 per 100k
    # per day, usually closer to 10 and then spiking exponentially in
    # December. This gives an average incidence slightly higher than the UK's.
    # I'm going to take 0.33 as a guess.
    # (LA Covid Data: https://www.nytimes.com/interactive/2021/us/los-angeles-california-covid-cases.html)
    source="https://www.rcgp.org.uk/getmedia/e564cc01-bc66-4165-aba5-c762b693f50d/2020-December.zip",
)


def estimate_prevalences():
    return [
        total_prevalence.__mul__(pandemic_decrease_factor),
    ]
