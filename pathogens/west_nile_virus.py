from pathogen_properties import *
from populations import us_population

# TODO: We should add a CI once we know how to get from tested individuals to
# CI

background = """“West Nile virus (WNV) is the leading cause of mosquito-borne 
disease in the continental United States.  It is most commonly spread to 
people by the bite of an infected mosquito. Cases of WNV occur during mosquito 
season, which starts in the summer and continues through fall. There are no 
vaccines to prevent or medications to treat WNV in people. About 1 in 5 people 
who are infected develop a fever and other symptoms” (CDC)."""


pathogen_chars = PathogenChars(
    na_type=NAType.RNA,
    enveloped=Enveloped.ENVELOPED,
    taxid=TaxID(11082),
)

LA_county_cases_in_2020 = IncidenceAbsolute(
    annual_infections=90,
    country="United States",
    state="California",
    county="Los Angeles County",
    date="2020",
    source="https://westnile.ca.gov/pdfs/VBDSAnnualReport20.pdf#?page=23",
)

# 3/4 of people with WNV are asymptomatic, and therefore probably did not get
# tested
asymptomatic_multiplier = Scalar(
    scalar=4,
    source="https://www.cdc.gov/mmwr/volumes/70/ss/ss7001a1.htm#:~:text=An%20estimated%2070%25%E2%80%9380%25%20of%20WNV%20infections%20are%20asymptomatic%20(8%2C9)",
)


def estimate_incidences() -> list[IncidenceRate]:
    return [
        LA_county_cases_in_2020.to_rate(
            us_population(
                county="Los Angeles County", state="California", year=2020
            )
        )
        * asymptomatic_multiplier
    ]


def estimate_prevalences() -> list[Prevalence]:
    return []
