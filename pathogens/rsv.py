from pathogen_properties import *
import math

background = """Respiratory syncytial virus, or RSV, is a common respiratory virus that usually causes mild, cold-like symptoms. Most people recover in a week or two, but RSV can be serious, especially for infants and older adults. RSV is the most common cause of bronchiolitis (inflammation of the small airways in the lung) and pneumonia (infection of the lungs) in children younger than 1 year of age in the United States."""

pathogen_chars = PathogenChars(
    na_type=NAType.RNA,
    enveloped=Enveloped.ENVELOPED,
    taxid=12814,
)

# positive tests in california during fall 2021
california_positive_tests_fall_2021 = Scalar(
        scalar = 6400000,
        country="United States",
        state = "California",
        date = "8/21 - 12/18, 2020",
        source="https://www.cdc.gov/surveillance/nrevss/images/rsvstate/RSV14NumCent5AVG_StateCA.htm",
)

# across how many days were the 6.4 million samples positive
testing_interval = Scalar(
        scalar = 119,
)    

# What percent of RSV cases remained during the pandemic, compared to pre/post-pandemic?
RSV_pandemic_decrease = Scalar(
        scalar = 0.14,
        state = "New York",
        source="https://www.dovepress.com/the-impact-of-the-covid-19-pandemic-on-respiratory-syncytial-virus-inf-peer-reviewed-fulltext-article-IDR#ref-cit0015",
)

# how many days does RSV usually last for
disease_duration = Scalar(
        scalar = 10,
        source = "https://www.cdc.gov/rsv/about/symptoms.html#:~:text=Most%20RSV%20infections%20go%20away%20on%20their%20own%20in%20a%20week%20or%20two.",
)    

CA_population = Scalar(
        scalar=40000000,
)    

def estimate_prevalences():
        return[
            RSV_pandemic_decrease.scalar * california_positive_tests_fall_2021.scalar * disease_duration.scalar / testing_interval.scalar / CA_population.scalar * 100000
         ]
