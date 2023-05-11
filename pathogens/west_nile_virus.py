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

west_nile_duration = SheddingDuration(
    # Symptoms last for 3-6 days usually, but sometimes for up to a month. I'm
    # going to use 7 days as an estimate, since I cannot find a better source
    # on this
    days=10,
    confidence_interval=(2, 18),
    source="https://myhealth.alberta.ca/Health/aftercareinformation/pages/conditions.aspx?hwid=abo5809#:~:text=In%20mild%20cases%20of%20West%20Nile%2C%20symptoms%20usually%20last%20for%203%20to%206%20days%2C%20and%20you%20can%20recover%20at%20home.%20If%20you%20get%20a%20more%20severe%20case%20of%20West%20Nile%2C%20symptoms%20can%20last%20for%20weeks%20or%20months%2C%20and%20you%20may%20need%20to%20stay%20in%20the%20hospital%20so%20you%20can%20get%20medicine%20to%20help%20you%20recover.",
)

# 3/4 of people with WNV are asymptomatic, and therefore probably did not get
# tested
asymptomatic_multiplier = Scalar(
    scalar=4,
    source="https://www.cdc.gov/mmwr/volumes/70/ss/ss7001a1.htm#:~:text=An%20estimated%2070%25%E2%80%9380%25%20of%20WNV%20infections%20are%20asymptomatic%20(8%2C9)",
)


def estimate_prevalences():
    return [
        LA_county_cases_in_2020.to_rate(
            us_population(
                county="Los Angeles County", state="California", year=2020
            )
        ).to_prevalence(west_nile_duration)
        * asymptomatic_multiplier
    ]
