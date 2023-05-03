from pathogen_properties import *

# ***We should add a CI once we know how to get from tested individuals to CI

background = """“West Nile virus (WNV) is the leading cause of mosquito-borne 
disease in the continental United States.  It is most commonly spread to 
people by the bite of an infected mosquito. Cases of WNV occur during mosquito 
season, which starts in the summer and continues through fall. There are no 
vaccines to prevent or medications to treat WNV in people. Fortunately, most 
people infected with WNV do not feel sick. About 1 in 5 people who are 
infected develop a fever and other symptoms” (CDC)."""


pathogen_chars = PathogenChars(
    na_type=NAType.RNA,
    enveloped=Enveloped.ENVELOPED,
    taxid=TaxID(11082),
)

LA_county_cases_in_2020 = IncidenceAbsolute(
    annual_infections=90,
    tag="LA-2020",
    source="https://westnile.ca.gov/pdfs/VBDSAnnualReport20.pdf#?page=23",
)

LA_county_population = Population(
    people=10_014_009,
    tag="LA-2020",
    source="https://www.census.gov/quickfacts/fact/table/losangelescountycalifornia,CA/POP010220#POP010220",
)


west_nile_duration = SheddingDuration(
    days=7,
    source="https://myhealth.alberta.ca/Health/aftercareinformation/pages/conditions.aspx?hwid=abo5809#:~:text=In%20mild%20cases%20of%20West%20Nile%2C%20symptoms%20usually%20last%20for%203%20to%206%20days%2C%20and%20you%20can%20recover%20at%20home.%20If%20you%20get%20a%20more%20severe%20case%20of%20West%20Nile%2C%20symptoms%20can%20last%20for%20weeks%20or%20months%2C%20and%20you%20may%20need%20to%20stay%20in%20the%20hospital%20so%20you%20can%20get%20medicine%20to%20help%20you%20recover.",
)

# 4/5 of people with WNV are asymptomatic, and therefore probably did not get
# tested
asymptomatic_multiplier = Scalar(
    scalar=5,
    source="https://www.cdc.gov/westnile/symptoms/index.html",
)


def estimate_prevalences():
    return [
        LA_county_cases_in_2020.to_rate(LA_county_population)
        .to_prevalence(west_nile_duration)
        .scale(asymptomatic_multiplier)
    ]
