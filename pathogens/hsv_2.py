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
    enveloped=Enveloped.ENVELOPED,
    taxid=TaxID(10310),
)


cdc_2018_nhanes_estimate = PrevalenceAbsolute(
    infections=18.6 * 1e6,
    confidence_interval=(18.1 * 1e6, 19.0 * 1e6),
    coverage_probability=0.5,
    date="2018",
    country="United States",
    tag="18-49yo",
    active=Active.LATENT,
    source="https://journals.lww.com/stdjournal/Fulltext/2021/04000/Estimates_of_the_Prevalence_and_Incidence_of.9.aspx#:~:text=In%202018%2C%20there,that%20are%20genital.",
    methods="https://journals.lww.com/stdjournal/Fulltext/2021/04000/Estimates_of_the_Prevalence_and_Incidence_of.9.aspx#:~:text=genital%20infections%20burden.-,METHODS,-We%20estimated%20the",
)

cdc_2015_2016_nhanes_estimate = Prevalence(
    infections_per_100k=0.121 * 100_000,
    country="United States",
    confidence_interval=(0.0966 * 100_000, 0.1495 * 100_000),
    # This CDC estimate is based on the following NHANES data:
    # https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/HSV_I.htm#:~:text=3710-,LBXHE2%20%2D%20Herpes%20Simplex%20Virus%20Type%202,-Variable%20Name%3A
    start_date="2015",
    end_date="2016",
    active=Active.LATENT,
    source="https://www.cdc.gov/nchs/data/databriefs/db304_table.pdf?#page=3",
    methods="https://www.cdc.gov/nchs/products/databriefs/db304.htm#:~:text=Estimates%20were%20calculated,p%20%3C%200.05.",
)


stratified_us_pop = {
    "15-19": Population(
        people=21_445_493,
        date="2018",
        country="United States",
        source="https://data.census.gov/table?q=2018+us+population&t=Civilian+Population",
        tag="15-19yo",
    ),
    "20-24": Population(
        people=21_717_962,
        date="2018",
        country="United States",
        source="https://data.census.gov/table?q=2018+us+population&t=Civilian+Population",
        tag="20-24yo",
    ),
    "25-34": Population(
        people=45_344_674,
        date="2018",
        country="United States",
        source="https://data.census.gov/table?q=2018+us+population&t=Civilian+Population",
        tag="25-34yo",
    ),
    "35-44": Population(
        people=41_498_453,
        date="2018",
        country="United States",
        source="https://data.census.gov/table?q=2018+us+population&t=Civilian+Population",
        tag="35-44yo",
    ),
    "45-54": Population(
        people=41_605_244,
        date="2018",
        country="United States",
        source="https://data.census.gov/table?q=2018+us+population&t=Civilian+Population",
        tag="45-54yo",
    ),
}

us_population_2018_18_to_49yo = Population(
    people=stratified_us_pop["15-19"].people * 0.6
    + stratified_us_pop["20-24"].people
    + stratified_us_pop["25-34"].people
    + stratified_us_pop["35-44"].people
    + stratified_us_pop["45-54"].people * 0.5,
    date="2018",
    country="United States",
    tag="18-49yo",
    inputs=list(stratified_us_pop.values()),
)


def estimate_prevalences() -> list[Prevalence]:
    us_2018 = cdc_2018_nhanes_estimate.to_rate(us_population_2018_18_to_49yo)

    # HSV_1 prevalence should be close to constant, so extrapolate from
    # 2018 to 2020 and 2021.
    return [
        dataclasses.replace(us_2018, date_source=Variable(date="2020")),
        dataclasses.replace(us_2018, date_source=Variable(date="2021")),
    ]


def estimate_incidences() -> list[IncidenceRate]:
    return []
