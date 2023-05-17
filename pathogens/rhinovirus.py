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
# Ideally, we would find the under 19 population instead to match our under 19
# prevalence estimate
la_county_under_18_population = Population(
    # 0.211 = the proportion of LA County residents under 18
    people=0.211
    * us_population(
        state="California", county="Los Angeles County", year=2020
    ).people,
    date="2020",
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
rhinovirus_0_4 = IncidenceRate(
    annual_infections_per_100k=113.2 * 100,
    # Table 2, Rhinoviruses.
    # "Annual isolation rates of respiratory viruses, Tecumseh Michigan,
    # USA, 1976-81. Actual rates per 1000 person year"
    source="https://doi.org/10.1017/S0950268800050779",
    start_date="1976",
    end_date="1981",
    country="United States",
    state="Michigan",
    county="Lenawee County",
)
rhinovirus_5_19 = IncidenceRate(
    annual_infections_per_100k=25.3 * 100,
    # Table 2, Rhinoviruses
    # "Annual isolation rates of respiratory viruses, Tecumseh Michigan,
    # USA, 1976-81. Actual rates per 1000 person year"
    source="https://doi.org/10.1017/S0950268800050779",
    start_date="1976",
    end_date="1981",
    country="United States",
    state="Michigan",
    county="Lenawee County",
)
rhinovirus_20_39 = IncidenceRate(
    annual_infections_per_100k=38.7 * 100,
    # Table 2, Rhinoviruses
    # "Annual isolation rates of respiratory viruses, Tecumseh Michigan,
    # USA, 1976-81. Actual rates per 1000 person year"
    source="https://doi.org/10.1017/S0950268800050779",
    start_date="1976",
    end_date="1981",
    country="United States",
    state="Michigan",
    county="Lenawee County",
)
rhinovirus_40_plus = IncidenceRate(
    annual_infections_per_100k=9.7 * 100,
    # Table 2, Rhinoviruses
    # "Annual isolation rates of respiratory viruses, Tecumseh Michigan,
    # USA, 1976-81. Actual rates per 1000 person year"
    source="https://doi.org/10.1017/S0950268800050779",
    start_date="1976",
    end_date="1981",
    country="United States",
    state="Michigan",
    county="Lenawee County",
)

la_county_0_4 = Population(
    people=285_140 + 273_131,
    date="2020",
    country="United States",
    state="California",
    county="Los Angeles County",
    source="http://proximityone.com/chartgraphics/pp06037_2020_001.htm",
)
la_county_5_19 = Population(
    people=(306_835 + 291_053) + (320_666 + 303_381) + (317_657 + 306_849),
    date="2020",
    country="United States",
    state="California",
    county="Los Angeles County",
    source="http://proximityone.com/chartgraphics/pp06037_2020_001.htm",
)
la_county_20_39 = Population(
    people=(330_183 + 327_625)
    + (406_008 + 399_410)
    + (413_566 + 394_103)
    + (370_948 + 355_002),
    date="2020",
    country="United States",
    state="California",
    county="Los Angeles County",
    source="http://proximityone.com/chartgraphics/pp06037_2020_001.htm",
)
la_county = Population(
    people=4_965_022 + 5_048_987,
    date="2020",
    country="United States",
    state="California",
    county="Los Angeles County",
    source="http://proximityone.com/chartgraphics/pp06037_2020_001.htm",
)

la_county_40_plus = (
    la_county - la_county_0_4 - la_county_5_19 - la_county_20_39
)

rhinovirus_old_yearround_study_estimate = (
    IncidenceRate.weightedAverageByPopulation(
        (rhinovirus_0_4, la_county_0_4),
        (rhinovirus_5_19, la_county_5_19),
        (rhinovirus_20_39, la_county_20_39),
        (rhinovirus_40_plus, la_county_40_plus),
    )
)

average_colds_per_year_adults = Scalar(
    scalar=2.5,
    country="United States",
    state="California",
    county="Los Angeles County",
    date="2023",
    # These numbers basically just seem reasonable to me and seem to be
    # generally reasonable according to the internet. I think this is more
    # accurate than trusting a very-old study with questionable methods
    source="https://my.clevelandclinic.org/health/diseases/12342-common-cold#:~:text=Adults%20catch%20two%20to%20three%20colds%20a%20year%2C%20while%20young%20children%20come%20down%20with%20a%20cold%20four%20or%20more%20times%20a%20year.",
)

average_colds_per_year_children = Scalar(
    scalar=4,
    country="United States",
    state="California",
    county="Los Angeles County",
    date="2023",
    # These numbers basically just seem reasonable to me and seem to be
    # generally reasonable according to the internet. I think this is more
    # accurate than trusting a very-old study with questionable methods
    source="https://my.clevelandclinic.org/health/diseases/12342-common-cold#:~:text=Adults%20catch%20two%20to%20three%20colds%20a%20year%2C%20while%20young%20children%20come%20down%20with%20a%20cold%20four%20or%20more%20times%20a%20year.",
)

annual_colds_la_county = IncidenceAbsolute(
    # "Adults catch two to three colds a year, while young children come down
    # with a cold four or more times a year."
    # Since we're concerned with fall months, I'm multiplying by 1.5x
    annual_infections=1.5
    * average_colds_per_year_adults.scalar
    * la_county_adult_population.people
    + 1.5
    * average_colds_per_year_children.scalar
    * la_county_under_18_population.people,
    country="United States",
    state="California",
    county="Los Angeles County",
    date="2020",
    # These numbers basically just seem reasonable to me and seem to be
    # generally reasonable according to the internet. I think this is more
    # accurate than trusting a very-old study with questionable methods
    source="https://my.clevelandclinic.org/health/diseases/12342-common-cold#:~:text=Adults%20catch%20two%20to%20three%20colds%20a%20year%2C%20while%20young%20children%20come%20down%20with%20a%20cold%20four%20or%20more%20times%20a%20year.",
)


annual_colds_per_100k = annual_colds_la_county.to_rate(
    us_population(state="California", county="Los Angeles County", year=2020)
)

# In the fall, rhinovirus accounts for a higher percentage of colds than
# during other times of year, according to this study
fall_proportion_of_colds_caused_by_rhinovirus = Scalar(
    scalar=0.82,
    source="https://journals.asm.org/doi/epdf/10.1128/jcm.35.11.2864-2868.1997?src=getftr",
)

fall_proportion_of_colds_caused_by_rhinovirus_estimate_2 = Scalar(
    scalar=(0.77 + 0.65 + 0.82) / 3,
    # This is a literature review from 2001 which gives a few different
    # estimates for what share of colds are from rhinoviruses during the fall
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7133757/",
)

overall_fall_proportion_of_colds_caused_by_rhinovirus = Scalar(
    # Average of the first estimate & 3 studies in the second estimate
    scalar=(0.82 + 0.77 + 0.65 + 0.82)
    / 4,
)


rhinovirus_shedding_duration = SheddingDuration(
    days=12,
    confidence_interval=(10, 14),
    source="https://erj.ersjournals.com/content/44/1/169#:~:text=Virus%20shedding%20lasts%20on%20average%20for%2010%E2%80%9314%20days%20in%20immunocompetent%20subjects",
)

rhinovirus_prevalence_using_colds = (
    annual_colds_per_100k.to_prevalence(rhinovirus_shedding_duration)
    * overall_fall_proportion_of_colds_caused_by_rhinovirus
)


finland_pandemic_decrease_factor = Scalar(
    scalar=1,
    source="https://onlinelibrary.wiley.com/doi/full/10.1002/jmv.27857#:~:text=Rhinovirus%20detections%20remained%20practically%20unchanged%20in%20all%20age%20groups%20throughout%20the%20pandemic%20period%20in%20Finland%20(Figure%C2%A02D).",
)

# Link for all RCGP reports (specific reference is given under source):
#  https://www.rcgp.org.uk/clinical-and-research/our-programmes/research-and-surveillance-centre/public-health-data
pandemic_decrease_factor = Scalar(
    scalar=0.4,
    date="2020",
    # See page 8 of the December 2020 report, around weeks 44 to 53
    # The 5-year avg at this time is in blue, and the national avg in red below
    # From this graph, the UK's decrease factor was around 0.2.
    # On page 5, we see that Covid incidence in the relevant weeks
    # is around 14-28 per 100k per day.
    # During fall 2020, LA County had an incidence of around 10-150 per 100k
    # per day, usually closer to 10 and then spiking exponentially in
    # December. This gives an average incidence slightly higher than the UK's.
    # Additionally, in Finland, Rhinovirus barely decreased at all.
    # I'm going to take 0.4 as a guess.
    # (LA Covid Data: https://www.nytimes.com/interactive/2021/us/los-angeles-california-covid-cases.html)
    source="https://www.rcgp.org.uk/getmedia/e564cc01-bc66-4165-aba5-c762b693f50d/2020-December.zip",
)


def estimate_prevalences():
    return [
        rhinovirus_old_yearround_study_estimate.to_prevalence(
            rhinovirus_shedding_duration
        )
        * (pandemic_decrease_factor),
        rhinovirus_prevalence_using_colds * (pandemic_decrease_factor),
    ]
