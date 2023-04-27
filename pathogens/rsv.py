from pathogen_properties import *

background = """Respiratory syncytial virus, or RSV, is a common respiratory virus that usually causes mild, cold-like symptoms.
Most people recover in a week or two, but RSV can be serious, especially for infants and older adults. 
RSV is the most common cause of bronchiolitis (inflammation of the small airways in the lung) and pneumonia (infection of the lungs)
 in children younger than 1 year of age in the United States."""

pathogen_chars = PathogenChars(
    na_type=NAType.RNA,
    enveloped=Enveloped.ENVELOPED,
    taxid=12814,
)

# tests were collected across 119 days, so multiply by 365/119 for yearly
california_positive_tests_fall_2021 = IncidenceAbsolute(
    annual_infections=6400000 * 365 / 119,
    country="United States",
    state="California",
    date="2020",
    start_date="2020-08-21",
    end_date="2020-12-18",
    source="https://www.cdc.gov/surveillance/nrevss/images/rsvstate/RSV14NumCent5AVG_StateCA.htm",
    tag="california fall 2020",
)

# What percent of RSV cases remained during the pandemic, compared to pre/post-pandemic?
RSV_pandemic_decrease = Scalar(
    scalar=0.14,
    state="New York",
    source="https://www.dovepress.com/the-impact-of-the-covid-19-pandemic-on-respiratory-syncytial-virus-inf-peer-reviewed-fulltext-article-IDR#ref-cit0015",
)

shedding_duration = SheddingDuration(
    days=14.1,
    source="https://academic.oup.com/aje/article/190/12/2536/6313422",
)

CA_population = Population(
    people=39538245,
    source="https://www.census.gov/quickfacts/CA",
    tag="california fall 2020",
)


def estimate_prevalences():
    return [
        california_positive_tests_fall_2021.to_rate(CA_population)
        .to_prevalence(shedding_duration)
        .scale(RSV_pandemic_decrease)
    ]
