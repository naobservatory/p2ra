from pathogen_properties import *

background = """Herpes Simplex Virus 2 is a common herpesvirus that causes
genital herpes and is largely sexually transmited. Most HSV-2 infections
persist lifelong (CDC, https://www.hopkinsmedicine.org/health/
conditions-and-diseases/herpes-hsv1-and-hsv2), manifesting as cold sores or
fever blisters (https://www.hopkinsmedicine.org/health/conditions-and-diseases/
herpes-hsv1-and-hsv2). After initial infection and potential symptoms, most
HSV-1 infections persist lifelong."""


pathogen_chars = PathogenChars(
    na_type=NAType.DNA,
    enveloped=Enveloped.NON_ENVELOPED,
    taxid=10310,
)


cdc_2015_2016_nhanes_seroprevalence = Prevalence(
    infections_per_100k=0.102 * 100_000,
    number_of_participants=3710,
    country="United States",
    start_date="2015",
    end_date="2016",
    source="https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/HSV_I.htm#:~:text=3710-,LBXHE2%20%2D%20Herpes%20Simplex%20Virus%20Type%202,-Variable%20Name%3A,",
)


cdc_2018_nhanes_estimate = PrevalenceAbsolute(
    infections=18.6 * 1e6,
    confidence_interval=(18.1 * 1e6, 19.0 * 1e6),
    date="2018",
    country="United States",
    source="https://journals.lww.com/stdjournal/Fulltext/2021/04000/Estimates_of_the_Prevalence_and_Incidence_of.9.aspx#:~:text=In%202018%2C%20there,that%20are%20genital.",
    methods="https://journals.lww.com/stdjournal/Fulltext/2021/04000/Estimates_of_the_Prevalence_and_Incidence_of.9.aspx#:~:text=genital%20infections%20burden.-,METHODS,-We%20estimated%20the",
)


cdc_2015_2016_nhanes_estimate = Prevalence(
    infections_per_100k=0.121 * 100_000,
    country="United States",
    confidence_interval=(0.0966 * 100_000, 0.1495 * 100_000),
    start_date="2015",
    end_date="2016",
    source="https://www.cdc.gov/nchs/data/databriefs/db304_table.pdf?#page=3",
    methods="https://www.cdc.gov/nchs/products/databriefs/db304.htm#:~:text=Estimates%20were%20calculated,p%20%3C%200.05.",
)


us_population_2018 = Population(
    people=327.2 * 1e6,
    country="United States",
    source="https://data.census.gov/table?q=2018+us+population&t=Civilian+Population",
    date="2018",
)


def estimate_prevalences():
    return [
        cdc_2015_2016_nhanes_seroprevalence,
        cdc_2015_2016_nhanes_estimate,
        cdc_2018_nhanes_estimate.to_rate(us_population_2018),
    ]
