from pathogen_properties import *
from populations import us_population

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

RHINOVIRUS_A = TaxID(147711)
RHINOVIRUS_B = TaxID(147712)
RHINOVIRUS_C = TaxID(463676)

pathogen_chars = PathogenChars(
    na_type=NAType.RNA,
    enveloped=Enveloped.NON_ENVELOPED,
    taxids=frozenset([RHINOVIRUS_A, RHINOVIRUS_B, RHINOVIRUS_C]),
)


la_county_under_18_population = Population(
    # 0.211 = the proportion of LA County residents under 18
    people=0.211
    * us_population(
        state="California", county="Los Angeles County", year=2020
    ).people,
    date="2020",
    country="United States",
    state="California",
    county="Los Angeles County",
    tag="under 18",
    source="https://www.census.gov/quickfacts/fact/table/losangelescountycalifornia#:~:text=Persons%20under-,18,-years%2C%20percent",
)

la_county_adult_population = Population(
    # 0.789 = the proportion of LA County residents over 18
    people=0.789
    * us_population(
        state="California", county="Los Angeles County", year=2020
    ).people,
    date="2020",
    country="United States",
    state="California",
    county="Los Angeles County",
    tag="over 18",
    source="https://www.census.gov/quickfacts/fact/table/losangelescountycalifornia#:~:text=Persons%20under-,18,-years%2C%20percent",
)


"""During this study, "A single family respondent was contacted weekly by 
telephone to obtain information on the onset of acute respiratory or enteric 
illnesses in any of the household members. Across this time they tracked the 
number of respiratory illnesses, across age groups. These are overall 
respiratory illnesses (which we might equate with colds). In the same study, 
they report that colds were isolated in 34% of cases. (note that across 13140 
cases of respiratory illnesses, they only sampled 2227 of cases for isolation, 
which could lead to bias in this isolation rate). But for now, let's use this
number. Checking a table of isolation rates across ages, the rate of
Rhinoviruses among infections roughly holds true across ages."""
annual_infections_per_100k_by_age_group = {
    "0-4": 113.2 * 100,
    "5-19": 25.3 * 100,
    "20-39": 38.7 * 100,
    "40+": 9.7 * 100,
}
rhinovirus_incidence_rates_by_age = {
    age_group: IncidenceRate(
        annual_infections_per_100k=incidence,
        start_date="1976",
        end_date="1981",
        country="United States",
        state="Michigan",
        county="Lenawee County",
        source="https://doi.org/10.1017/S0950268800050779",
    )
    for age_group, incidence in annual_infections_per_100k_by_age_group.items()
}
people_per_age_group = {
    "0-4": 285_140 + 273_131,
    "5-19": (306_835 + 291_053) + (320_666 + 303_381) + (317_657 + 306_849),
    "20-39": (330_183 + 327_625)
    + (406_008 + 399_410)
    + (413_566 + 394_103)
    + (370_948 + 355_002),
    "overall": 4_965_022 + 5_048_987,
}

people_per_age_group["40+"] = (
    people_per_age_group["overall"]
    - people_per_age_group["0-4"]
    - people_per_age_group["5-19"]
    - people_per_age_group["20-39"]
)

la_county_populations_by_age = {
    age_group: Population(
        people=numPeople,
        date="2020",
        country="United States",
        state="California",
        county="Los Angeles County",
        source="http://proximityone.com/chartgraphics/pp06037_2020_001.htm",
    )
    for age_group, numPeople in people_per_age_group.items()
}

rhinovirus_1970s_tecumseh_study_la_estimate = (
    IncidenceRate.weightedAverageByPopulation(
        (
            rhinovirus_incidence_rates_by_age["0-4"],
            la_county_populations_by_age["0-4"],
        ),
        (
            rhinovirus_incidence_rates_by_age["5-19"],
            la_county_populations_by_age["5-19"],
        ),
        (
            rhinovirus_incidence_rates_by_age["20-39"],
            la_county_populations_by_age["20-39"],
        ),
        (
            rhinovirus_incidence_rates_by_age["40+"],
            la_county_populations_by_age["40+"],
        ),
    )
)

average_colds_per_year_adults = IncidenceRate(
    annual_infections_per_100k=2.5 * 100_000,
    country="United States",
    date="2023",
    # "Adults catch two to three colds a year, while young children come down
    # with a cold four or more times a year."
    # These numbers basically just seem reasonable to me and seem to be
    # generally reasonable according to the internet. I think this is more
    # accurate than trusting a very-old study with questionable methods
    source="https://my.clevelandclinic.org/health/diseases/12342-common-cold#:~:text=Adults%20catch%20two%20to%20three%20colds%20a%20year%2C%20while%20young%20children%20come%20down%20with%20a%20cold%20four%20or%20more%20times%20a%20year.",
)

average_colds_per_year_children = IncidenceRate(
    annual_infections_per_100k=4 * 100_000,
    country="United States",
    date="2023",
    # "Adults catch two to three colds a year, while young children come down
    # with a cold four or more times a year."
    # These numbers basically just seem reasonable to me and seem to be
    # generally reasonable according to the internet. I think this is more
    # accurate than trusting a very-old study with questionable methods
    source="https://my.clevelandclinic.org/health/diseases/12342-common-cold#:~:text=Adults%20catch%20two%20to%20three%20colds%20a%20year%2C%20while%20young%20children%20come%20down%20with%20a%20cold%20four%20or%20more%20times%20a%20year.",
)


colds_la_county = IncidenceRate.weightedAverageByPopulation(
    (average_colds_per_year_children, la_county_under_18_population),
    (average_colds_per_year_adults, la_county_adult_population),
)


# In the fall, rhinovirus accounts for a higher percentage of colds than
# during other times of year
fall_proportion_of_colds_caused_by_rhinovirus = Scalar(
    scalar=(0.77 + 0.65 + 0.82) / 3,
    # This is a literature review from 2001 which gives a few different
    # estimates for what share of colds are from rhinoviruses during the fall
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7133757/",
)


# Rhinovirus incidence did not significantly change during the pandemic, so I
# think we can apply these post-pandemic numbers to fall 2020 without further
# adjustment.
# hRV did not decrease much in Canada: https://onlinelibrary.wiley.com/doi/10.1111/irv.12930#:~:text=The%20causes%20of,when%20PHMs%20ease.
# hRV also did not decrease much in Guangzhou, China:
# https://onlinelibrary.wiley.com/doi/10.1002/iid3.381#:~:text=A%20total%20of%2061%20samples%20were%20collected%20and%2017%20(27.87%25)%20of%20which%20were%20positive%20for%20rhinovirus
# This study finds that hRV actually increased slightly in Australia: https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.47.2001847#:~:text=In%20contrast%2C%20rhinovirus%20detections%20were%20well%20above%20average.
# hRV also increased slightly in Japan: https://onlinelibrary.wiley.com/doi/10.1111/irv.12854#:~:text=the%20frequency%20of%20rhinovirus%20infection%20increased%20appreciably%20during%20the%20COVID%2D19%20pandemic
# hRV also did not decrease in Finland - This source could also be used
# as a source for rhinovirus incidence: "https://onlinelibrary.wiley.com/doi/full/10.1002/jmv.27857#:~:text=Rhinovirus%20detections%20remained%20practically%20unchanged%20in%20all%20age%20groups%20throughout%20the%20pandemic%20period%20in%20Finland%20(Figure%C2%A02D)."


def estimate_incidences():
    return (
        rhinovirus_1970s_tecumseh_study_la_estimate,
        colds_la_county * fall_proportion_of_colds_caused_by_rhinovirus,
    )


def estimate_prevalences():
    return ()
