from pathogen_properties import *

background = """Hepatitis B is a liver infection caused by the Hepatitis B
 virus. It is transmitted through birth and contact with infected blood or
 bodily fluids. Hepatitis B infection increases the risk for hepatocellular
 carcinoma"""

# TODO: Add in-house NHANES estimate. Specifically, we could present NHANES data
# for HBV-core antibody (evidence of past infection), and HBV-surface antigen
# (evidence of active acute or chronic infection):
# https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/HEPBD_J.htm


pathogen_chars = PathogenChars(
    na_type=NAType.DNA,
    enveloped=Enveloped.ENVELOPED,
    taxid=TaxID(10407),
)


hbv_us_2022_deaths_wonder = Number(
    # Would need to find CFR to arrive at prevalence
    number=1_240,
    date="2022",
    country="United States",
    source="https://wonder.cdc.gov/controller/saved/D176/D341F095",
)

dna_present_in_serum = SheddingDuration(
    days=2 * 30.5,  # 2 months
    date="1991",
    source="https://pubmed.ncbi.nlm.nih.gov/1909458/",  # paywalled, data for
    # this variable was read from figure 2, showing duration of presence of
    # HBV-DNA in blood.
)

cdc_estimated_acute_2019 = IncidenceAbsolute(
    annual_infections=20_700,
    confidence_interval=(11_800, 50_800),  # 95% Bootstrap Confidence Interval
    coverage_probability=0.95,
    country="United States",
    date="2019",
    tag="us_2019",
    source="https://www.cdc.gov/hepatitis/statistics/2019surveillance/Introduction.htm#Technical:~:text=20%2C700%20estimated%20infections%20(95%25%20CI%3A%2011%2C800%E2%80%9350%2C800)",
    methods="https://www.cdc.gov/hepatitis/statistics/2019surveillance/Introduction.htm#Technical:~:text=To%20account%20for,CI%3A%2011.0%E2%80%9347.4).",
)

estimated_chronic_us_2020 = PrevalenceAbsolute(
    # This estimate is higher than what NHANES data alone would suggest. This
    # is because the study takes into account the high burden of HBV in
    # foreign-born individuals, and high-risk populations, both of which are
    # likely underrepresented in NHANES.
    infections=1_590_000,
    confidence_interval=(1_250_000, 2_490_000),
    active=Active.LATENT,
    country="United States",
    date="2020",
    tag="us_2020",
    source="https://journals.lww.com/ajg/fulltext/2020/09000/prevalence_of_chronic_hepatitis_b_virus_infection.20.aspx",
    # Methods for prevalence aggregation isn't given:
    methods="https://journals.lww.com/ajg/fulltext/2020/09000/prevalence_of_chronic_hepatitis_b_virus_infection.20.aspx#:~:text=Panel%20members%20researched,in%20the%20US.",
)


us_population_2019 = Population(
    people=328.2 * 1e6,
    country="United States",
    date="2019",
    tag="us_2019",
    source="https://data.census.gov/table?q=us+population+2019&tid=ACSDP1Y2019.DP05",
)


us_population_2020 = Population(
    people=331.4 * 1e6,
    country="United States",
    date="2020",
    tag="us_2020",
    source="https://data.census.gov/table?q=2020+us+population&t=Civilian+Population&tid=DECENNIALPL2020.P1",
)


def estimate_prevalences():
    return [
        cdc_estimated_acute_2019.to_rate(us_population_2019).to_prevalence(
            dna_present_in_serum
        ),
        estimated_chronic_us_2020.to_rate(us_population_2020),
    ]
