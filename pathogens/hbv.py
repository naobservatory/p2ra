from pathogen_properties import *
from populations import us_population

background = """Hepatitis B is a liver infection caused by the Hepatitis B
 virus. It is transmitted through birth and contact with infected blood or
 bodily fluids. Hepatitis B infection increases the risk for hepatocellular
 carcinoma"""

# TODO:
# - Add in-house NHANES estimate. Specifically, we could present NHANES data
# for HBV-core antibody (evidence of past infection), and HBV-surface antigen
# (evidence of active acute or chronic infection):
# https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/HEPBD_J.htm
# - Incorporate an estimate based on mortality date. This data can be acquired
# from the CDC WONDER database (https://wonder.cdc.gov/mcd-icd10-provisional.html). Parameters for relevant query is:
# - Select time period of death": 2022, and
# - Slect underlying cause of death: B16 (Acute hepatitis B).
# All other fields can be left as default.


pathogen_chars = PathogenChars(
    na_type=NAType.DNA,
    enveloped=Enveloped.ENVELOPED,
    taxid=TaxID(10407),
)

cdc_estimated_acute_2019 = IncidenceAbsolute(
    annual_infections=20_700,
    # During 2019, a total of 3,192 acute hepatitis B cases were reported to
    # CDC, resulting in 20,700 estimated infections (95% CI: 11,800–50,800)
    # after adjusting for case underascertainment and underreporting
    confidence_interval=(11_800, 50_800),  # 95% Bootstrap Confidence Interval
    coverage_probability=0.95,
    country="United States",
    # I picked 2019 rates, as estimated rates in 2020 were 30% lower, even
    # while deaths stayed the same: https://www.cdc.gov/hepatitis/statistics/2020surveillance/introduction/national-profile.htm#:~:text=hepatitis%20B%20transmission.-,Data%20from%20death%20certificates%20filed%20in%20the%20vital%20records%20offices%20of,same%20as%20the%20rate%20during%202019%20(0.42%20deaths%20per%20100%2C000%20population).,-Hepatitis%20C
    date="2019",
    source="https://www.cdc.gov/hepatitis/statistics/2019surveillance/Introduction.htm#Technical:~:text=20%2C700%20estimated%20infections%20(95%25%20CI%3A%2011%2C800%E2%80%9350%2C800)",
    methods="https://www.cdc.gov/hepatitis/statistics/2019surveillance/Introduction.htm#Technical:~:text=To%20account%20for,CI%3A%2011.0%E2%80%9347.4).",
)

estimated_chronic_us_2020 = PrevalenceAbsolute(
    infections=1_721_027,
    # This is their mid estimate in table 4, second to last row. The estimate
    # listed in their paper and on the last line of the table doesn't add up.
    # They say they added 303'237 (point estimate US-born) + 1'076'069 (lower
    # estimate Foreign-born) + 72'013 (point estimate Non-NHANES). But that
    # adds up to 1'450'320, which isn't the number they give (1.59M). We are
    # contacting the authors on this.
    confidence_interval=(1_249_055, 2_491_435),  # These numbers represent a
    # "low" and "high" estimate, labels that aren't further defined. Numbers
    # can be found in table 4, second to last row.
    active=Active.LATENT,
    country="United States",
    date="2020",
    source="https://journals.lww.com/ajg/fulltext/2020/09000/prevalence_of_chronic_hepatitis_b_virus_infection.20.aspx#:~:text=Table%204.%3A%20Estimated%20prevalence%20of%20chronic%20HBV%20in%20the%20United%20Statesa",
    methods="https://journals.lww.com/ajg/fulltext/2020/09000/prevalence_of_chronic_hepatitis_b_virus_infection.20.aspx#:~:text=Panel%20members%20researched,in%20the%20US.",
)


us_population_2019 = Population(
    people=328.2 * 1e6,
    country="United States",
    date="2019",
    source="https://data.census.gov/table?q=us+population+2019&tid=ACSDP1Y2019.DP05",
)


ohio_acute_incidence_2021 = IncidenceRate(
    # This source reports both an acute and a chronic incidence. I think it
    # makes sense to add these 2 numbers for an acute estimate, since chronic
    # cases still have an acute phase in which they shed.
    annual_infections_per_100k=16.5,
    country="United States",
    state="Ohio",
    date="2021",
    # Note that this source has yearly data for 2017-2021 if it's ever needed
    source="https://odh.ohio.gov/know-our-programs/viral-hepatitis/data-statistics/hbv_5yr_report",
)


def estimate_prevalences() -> list[Prevalence]:
    return [estimated_chronic_us_2020.to_rate(us_population(year=2020))]


def estimate_incidences() -> list[IncidenceRate]:
    return [
        cdc_estimated_acute_2019.to_rate(us_population_2019),
        ohio_acute_incidence_2021,
    ]
