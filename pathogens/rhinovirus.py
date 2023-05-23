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

# FWIW, this study claims that Covid Incidence has a significant positive
# correlation with rhinovirus incidence: https://onlinelibrary.wiley.com/doi/10.1111/irv.12930#:~:text=Figure%C2%A01%20suggests%20that%2C%20with%20the%20hindsight%20of%2018%E2%80%89months%20of%20observations%20for%20both%20hRV/EV%20and%20SARS%2DCoV%2D2%20infections%20in%20Canada%2C%20hRV/EV%20could%20have%20been%20used%20as%20a%20gauge%20for%20PHMs%20effectiveness%20as%20well%20as%20early%20warning%20for%20SARS%2DCoV%2D2%20resurgences]

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


rhinovirus_shedding_duration = SheddingDuration(
    # This estimate is for nasal shedding
    days=12,
    confidence_interval=(10, 14),
    source="https://erj.ersjournals.com/content/44/1/169#:~:text=Virus%20shedding%20lasts%20on%20average%20for%2010%E2%80%9314%20days%20in%20immunocompetent%20subjects",
)


la_county_under_18_population = Population(
    # 0.211 = the proportion of LA County residents under 18
    people=0.211
    * us_population(state="California", county="Los Angeles County", year=2020).people,
    date="2020",
    tag="under 18",
    source="https://www.census.gov/quickfacts/fact/table/losangelescountycalifornia#:~:text=Persons%20under-,18,-years%2C%20percent",
)

la_county_adult_population = Population(
    # 0.789 = the proportion of LA County residents over 18
    people=0.789
    * us_population(state="California", county="Los Angeles County", year=2020).people,
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

la_county_40_plus = la_county - la_county_0_4 - la_county_5_19 - la_county_20_39

rhinovirus_old_yearround_study_estimate = IncidenceRate.weightedAverageByPopulation(
    (rhinovirus_0_4, la_county_0_4),
    (rhinovirus_5_19, la_county_5_19),
    (rhinovirus_20_39, la_county_20_39),
    (rhinovirus_40_plus, la_county_40_plus),
)

#  "Adults catch two to three colds a year, while young children come down
# with a cold four or more times a year."
# These numbers basically just seem reasonable to me and seem to be
# generally reasonable according to the internet. I think this is more
# accurate than trusting a very-old study with questionable methods
average_colds_per_year_adults = IncidenceRate(
    annual_infections_per_100k=2.5 * 100_000,
    country="United States",
    date="2023",
    # These numbers basically just seem reasonable to me and seem to be
    # generally reasonable according to the internet. I think this is more
    # accurate than trusting a very-old study with questionable methods
    source="https://my.clevelandclinic.org/health/diseases/12342-common-cold#:~:text=Adults%20catch%20two%20to%20three%20colds%20a%20year%2C%20while%20young%20children%20come%20down%20with%20a%20cold%20four%20or%20more%20times%20a%20year.",
)

average_colds_per_year_children = IncidenceRate(
    annual_infections_per_100k=4 * 100_000,
    country="United States",
    date="2023",
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

# People get more colds in the fall. Since we are only interested in the fall,
# it makes sense to increase the yearround number by some reasonable factor
# Assuming people get no colds in the summer and an equal amount of colds
# during the fall, winter, and spring, this factor should be about 1.33
colds_fall_disproportionality = Scalar(scalar=1.33)

fall_rhinovirus_incidence_la_county = (
    colds_la_county * colds_fall_disproportionality
) * fall_proportion_of_colds_caused_by_rhinovirus

rhinovirus_prevalence_using_colds = fall_rhinovirus_incidence_la_county.to_prevalence(
    rhinovirus_shedding_duration
)


# hRV did not decrease much in Canada: https://onlinelibrary.wiley.com/doi/10.1111/irv.12930#:~:text=The%20causes%20of,when%20PHMs%20ease.
# hRV also did not decrease much in Guangzhou, China:
# https://onlinelibrary.wiley.com/doi/10.1002/iid3.381#:~:text=A%20total%20of%2061%20samples%20were%20collected%20and%2017%20(27.87%25)%20of%20which%20were%20positive%20for%20rhinovirus
# This study finds the hRV actually increased in Australia: https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.47.2001847#:~:text=In%20contrast%2C%20rhinovirus%20detections%20were%20well%20above%20average.
# hRV also increased in Japan: https://onlinelibrary.wiley.com/doi/10.1111/irv.12854#:~:text=the%20frequency%20of%20rhinovirus%20infection%20increased%20appreciably%20during%20the%20COVID%2D19%20pandemic
finland_china_japan_australia_canada_pandemic_decrease_factor = Scalar(
    scalar=1,
    # This source finds rhinovirus did not decrease during the pandmeic in
    # Finland - Could also be used as a source for rhinovirus incidence
    source="https://onlinelibrary.wiley.com/doi/full/10.1002/jmv.27857#:~:text=Rhinovirus%20detections%20remained%20practically%20unchanged%20in%20all%20age%20groups%20throughout%20the%20pandemic%20period%20in%20Finland%20(Figure%C2%A02D).",
)

# Given the many studies above that find that rhinovirus did not decrease
# during the pandemic, I got rid of the pandemic_decrease_factor scalar


def estimate_prevalences():
    return [
        rhinovirus_old_yearround_study_estimate.to_prevalence(
            rhinovirus_shedding_duration
        ),
        rhinovirus_prevalence_using_colds,
    ]
