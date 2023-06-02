from pathogen_properties import *
from populations import us_population

background = """Rhinoviruses are the most frequent cause of the common cold.
The virus is not as seasonal as viruses like influenza or coronviruses.
(https://journals.plos.org/plosone/article/figure?id=10.1371/journal.pone.
0114871.g001)"""


# TODOs for SIMON:
#  - Check if national estimates can be applied to regions like Ohio during
#  some specific time period.
#  - Check if Tecumseh based estimate should be adjusted for seasonality,
#  similar to the colds_la_county estimate.

RHINOVIRUS_A = TaxID(147711)
RHINOVIRUS_B = TaxID(147712)
RHINOVIRUS_C = TaxID(463676)

pathogen_chars = PathogenChars(
    na_type=NAType.RNA,
    enveloped=Enveloped.NON_ENVELOPED,
    taxids=frozenset([RHINOVIRUS_A, RHINOVIRUS_B, RHINOVIRUS_C]),
)


la_fraction_u18 = Scalar(
    scalar=0.211,
    source="https://www.census.gov/quickfacts/fact/table/losangelescountycalifornia/AGE135221#AGE135221",
)
la_population_u18 = (
    us_population(year=2020, state="California", county="Los Angeles County")
    * la_fraction_u18
)

la_population_18plus = (
    us_population(year=2020, state="California", county="Los Angeles County")
    - la_population_u18
)


TECUMSEH_SURVEILLANCE_STUDY_1993 = "https://doi.org/10.1017/S0950268800050779"

# During this study, run from 1965 to 1971 and 1976 to 1981, a single family
# respondent was contacted weekly by telephone to obtain information on the onset
# of acute respiratory or enteric illnesses in any of the household members.
# Across this time they tracked the number of respiratory illnesses, across age
# groups.  In the same study, they report that Rhinovirus was isolated in 34% of
# cases. (note that across 13140 cases of respiratory illnesses, they only
# sampled 2227 of cases for isolation, which could lead to bias in this
# isolation rate).

tecumseh_rhinovirus_infections_by_age = {
    # % rhinovirus infections * # colds / person-years * 100'000
    "0-4": 0.39 * (2657 / 539) * 100_000,
    "5-19": 0.21 * (4373 / 1541) * 100_000,
    "20-39": 0.5 * (3326 / 1523) * 100_000,
    "40+": 0.36 * (2784 / 1757) * 100_000,
}

rhinovirus_incidence_rates_by_age = {}

for age_group, incidence in tecumseh_rhinovirus_infections_by_age.items():
    incidence_rate = IncidenceRate(
        annual_infections_per_100k=incidence,
        start_date="1976",
        end_date="1981",
        country="United States",
        state="Michigan",
        county="Lenawee County",
        source="TECUMSEH_SURVEILLANCE_STUDY_1993 ",
    )
    rhinovirus_incidence_rates_by_age[age_group] = incidence_rate

LA_DEMOGRAPHIC_DATA_2020 = (
    "http://proximityone.com/chartgraphics/pp06037_2020_001.htm"
)

la_age_groups = {
    "0-4": 285_140 + 273_131,
    "5-19": (306_835 + 291_053) + (320_666 + 303_381) + (317_657 + 306_849),
    "20-39": (330_183 + 327_625)
    + (406_008 + 399_410)
    + (413_566 + 394_103)
    + (370_948 + 355_002),
}


la_county_populations_by_age = {}

for age_group, n_people in la_age_groups.items():
    age_group_population = Population(
        people=n_people,
        date="2020",
        country="United States",
        state="California",
        county="Los Angeles County",
        source=LA_DEMOGRAPHIC_DATA_2020,
    )
    la_county_populations_by_age[age_group] = age_group_population


under_40_population = Population(
    people=sum(la_county_populations_by_age.values()),
    date="2020",
    country="United States",
    state="California",
    county="Los Angeles County",
    source=LA_DEMOGRAPHIC_DATA_2020,
)


la_county_populations_by_age["40+"] = (
    us_population(year=2020, state="California", county="Los Angeles County")
    - under_40_population
)


rhinovirus_1970s_tecumseh_based_la_estimate = (
    # This calculation rests on the assumption that people still get the same
    # number of annual rhinoviruses for their age as they did in Tecumseh in
    # the 70s. Here, we merely adjust for the the difference in age breakdown
    # between the two contexts.
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
    source="https://my.clevelandclinic.org/health/diseases/12342-common-cold#:~:text=Adults%20catch%20two%20to%20three%20colds%20a%20year%2C%20while%20young%20children%20come%20down%20with%20a%20cold%20four%20or%20more%20times%20a%20year.",
)

average_colds_per_year_children = IncidenceRate(
    annual_infections_per_100k=4 * 100_000,
    country="United States",
    date="2023",
    source="https://my.clevelandclinic.org/health/diseases/12342-common-cold#:~:text=Adults%20catch%20two%20to%20three%20colds%20a%20year%2C%20while%20young%20children%20come%20down%20with%20a%20cold%20four%20or%20more%20times%20a%20year.",
)


colds_la_county = IncidenceRate.weightedAverageByPopulation(
    (average_colds_per_year_children, la_population_u18),
    (average_colds_per_year_adults, la_population_18plus),
)


fall_proportion_of_colds_caused_by_rhinovirus = Scalar(
    # "Rhinoviruses comprise more than three quarters of viruses circulating in early
    # autumn."
    scalar=0.75,
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7133757/",
)


# Rhinovirus incidence did not significantly change during the pandemic.
# Study from Canada: https://onlinelibrary.wiley.com/doi/10.1111/irv.12930#:~:text=The%20causes%20of,when%20PHMs%20ease.
# Study from Guangzhou, China:
# https://onlinelibrary.wiley.com/doi/10.1002/iid3.381#:~:text=A%20total%20of%2061%20samples%20were%20collected%20and%2017%20(27.87%25)%20of%20which%20were%20positive%20for%20rhinovirus
# Study from Australia: https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.47.2001847#:~:text=In%20contrast%2C%20rhinovirus%20detections%20were%20well%20above%20average.
# Study from Japan: https://onlinelibrary.wiley.com/doi/10.1111/irv.12854#:~:text=the%20frequency%20of%20rhinovirus%20infection%20increased%20appreciably%20during%20the%20COVID%2D19%20pandemic
# Study from Finland: https://onlinelibrary.wiley.com/doi/full/10.1002/jmv.27857#:~:text=Rhinovirus%20detections%20remained%20practically%20unchanged%20in%20all%20age%20groups%20throughout%20the%20pandemic%20period%20in%20Finland%20(Figure%C2%A02D)
# The final study could also be used as a source for rhinovirus incidence


def estimate_incidences():
    return (
        rhinovirus_1970s_tecumseh_based_la_estimate,
        colds_la_county * fall_proportion_of_colds_caused_by_rhinovirus,
    )


def estimate_prevalences():
    return ()
