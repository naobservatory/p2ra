from pathogen_properties import *

background = """Herpes Simplex Virus 2 is a common herpesvirus that causes genital herpes and is largely sexually transmited. Most HSV-2 infections persist lifelong (CDC, https://www.
hopkinsmedicine.org/health/conditions-and-diseases/herpes-hsv1-and-hsv2), manifesting as cold sores or fever
blisters (https://www.hopkinsmedicine.org/health/conditions-and-diseases/herpes-hsv1-and-hsv2). After initial
infection and potential symptoms, most HSV-1 infections persist lifelong."""


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
cdc_2015_2016_nhanes_estimate = Prevalence(
    infections_per_100k=0.121 * 100_000,
    country="United States",
    confidence_interval=(0.0966 * 100_000, 0.1495 * 100_000),
    country="United States",
    start_date="2015",
    end_date="2016",
    source="https://www.cdc.gov/nchs/data/databriefs/db304_table.pdf?#page=3",
    methods="https://www.cdc.gov/nchs/products/databriefs/db304.htm#:~:text=Estimates%20were%20calculated,p%20%3C%200.05.",
)


def estimate_prevalences():
    return [cdc_2015_2016_nhanes_seroprevalence, cdc_2015_2016_nhanes_estimate]
