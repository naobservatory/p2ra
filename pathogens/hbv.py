from pathogen_properties import *

background = """Hepatitis B is a liver infection caused by the hepatitis B virus. It is transmitted through birth and contact with infected blood or bodily fluids. Hepatitis B infection increases the risk for hepatocellular carcinoma"""

# TODO: Add case surveillance data: https://wonder.cdc.gov/nndss/static/2019/annual/2019-table2h.html


pathogen_chars = PathogenChars(
    na_type=NAType.DNA,
    enveloped=Enveloped.ENVELOPED,
    taxid=TaxID(10407),
)

# To add: in-house NHANES estimate. Specifically, NHANES tests for HBV-core
# antibody (evidence of past vaccination), HBV-core antibody (evidence of past
# infection), and HBV-surface antigen (evidence of active acute or chronic
# infection). Importantly surface antigen only gets tested if core antibody
# testing returns positive.
# https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/HEPBD_J.htm

hbv_us_2022_deaths_wonder = Number(
    # Would need to find CFR to arrive at prevalence
    number=1_240,
    date="2022",
    country="United States",
    source="https://wonder.cdc.gov/controller/saved/D176/D341F095",
)

estimated_acute_2019 = PrevalenceAbsolute(
    infections=20_700,
    confidence_interval=(11_800, 50_800),
    # 95% Bootstrap Confidence Interval
    coverage_probability=0.95,
    country="United States",
    date="2019",
    tag="us_2019",
    active=Active.ACTIVE,
    # Still need to find the specific paper this estimate is based on
    source="https://www.cdc.gov/hepatitis/statistics/2019surveillance/HepB.htm",
)

estimated_chronic_us_2020 = PrevalenceAbsolute(
    # This estimate is higher than what NHANES data alone would suggest.This is
    # because the study takes into account the high burden of HBV in
    # foreign-born individuals, and high-risk populations, both of which are
    # liklely underrepresented in NHANES.
    infections=1_590_000,
    confidence_interval=(1_250_000, 2_490_000),
    active=Active.LATENT,
    country="United States",
    date="2020",
    tag="us_2020",
    source="https://journals.lww.com/ajg/fulltext/2020/09000/prevalence_of_chronic_hepatitis_b_virus_infection.20.aspx",
    # Note that this study looks like an informal literature review:
    methods="https://journals.lww.com/ajg/fulltext/2020/09000/prevalence_of_chronic_hepatitis_b_virus_infection.20.aspx#:~:text=Panel%20members%20researched,in%20the%20US.",
)


# US population 2019

us_population_2019 = Population(
    people=328.2 * 1e6,
    country="United States",
    date="2019",
    tag="us_2019",
    source="https://data.census.gov/table?q=us+population+2019&tid=ACSDP1Y2019.DP05",
)

# US population 2020

us_population_2020 = Population(
    people=331.4 * 1e6,
    country="United States",
    date="2020",
    tag="us_2020",
    source="https://data.census.gov/table?q=2020+us+population&t=Civilian+Population&tid=DECENNIALPL2020.P1",
)


def estimate_prevalences():
    return [
        estimated_acute_2019.to_rate(us_population_2019),
        estimated_chronic_us_2020.to_rate(us_population_2020),
    ]
