from pathogen_properties import *
from populations import us_population

background = """Hepatitis B is a liver infection caused by the Hepatitis B
 virus. It is transmitted through birth and contact with infected blood or
 bodily fluids. Hepatitis B infection increases the risk for hepatocellular
 carcinoma"""

# TODO:
# 1. Add in-house NHANES estimate. Specifically, we could present NHANES data
# for HBV-core antibody (evidence of past infection), and HBV-surface antigen
# (evidence of active acute or chronic infection):
# https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/HEPBD_J.htm
# 2. Incorporate an estimate based on mortality date. This data can be acquired
# from the CDC WONDER database (https://wonder.cdc.gov/mcd-icd10-provisional.html).
# Parameters for relevant query is:
# - Select time period of death": 2022, and
# - Select underlying cause of death: B16 (Acute hepatitis B).
# All other fields can be left as default.


pathogen_chars = PathogenChars(
    na_type=NAType.DNA,
    enveloped=Enveloped.ENVELOPED,
    taxid=TaxID(10407),
    selection=SelectionRound.ROUND_1,
)

denmark_estimated_chronic_2007 = Prevalence(
    infections_per_100k=0.0024 * 100_000,
    # "the total national estimate was 10,668 (95% CI: 10,224–16,164),
    # corresponding to a prevalence of 0.24% (95% CI: 0.23–0.37%) in the
    # Danish population older than 15 years."
    confidence_interval=(0.0023 * 100_000, 0.0037 * 100_000),  # 95% CI
    coverage_probability=0.95,
    active=Active.LATENT,
    country="Denmark",
    date="2007",
    source="https://www.eurosurveillance.org/content/10.2807/1560-7917.ES2013.18.47.20637#:~:text=the%20total%20national%20estimate%20was%2010%2C668%20(95%25%20CI%3A%2010%2C224%E2%80%9316%2C164)%2C%20corresponding%20to%20a%20prevalence%20of%200.24%25%20(95%25%20CI%3A%200.23%E2%80%930.37%25)%20in%20the%20Danish%20population%20older%20than%2015%20years.",
)

us_population_2019 = Population(
    people=328.2 * 1e6,
    country="United States",
    date="2019",
    source="https://data.census.gov/table?q=us+population+2019&tid=ACSDP1Y2019.DP05",
)

cdc_acute_underreporting_factor_2020 = Scalar(
    scalar=6.5,
    confidence_interval=(3.7, 15.9),
    coverage_probability=0.95,
    country="United States",
    source="https://www.cdc.gov/hepatitis/statistics/2020surveillance/introduction/technical-notes.htm#:~:text=each%20reported%20case%20of%20acute%20hepatitis%20B%20represents%206.5%20estimated%20infections%20(95%25%20bootstrap%20CI%3A%203.7%E2%80%9315.9",
)


ohio_reported_acute_2021 = IncidenceRate(
    # This source reports both acute and total (acute+chronic) reported
    # incidence. We only report acute incidence here, as chronic prevalence is
    # covered by the national estimate.
    annual_infections_per_100k=0.9,
    country="United States",
    state="Ohio",
    date="2021",
    # Note that this source has yearly data for 2017-2021 if it's ever needed
    source="https://odh.ohio.gov/know-our-programs/viral-hepatitis/data-statistics/hbv_5yr_report",
)

estimated_acute_2019 = IncidenceAbsolute(
    annual_infections=20_700,
    # During 2019, a total of 3,192 acute hepatitis B cases were reported to
    # CDC, resulting in 20,700 estimated infections (95% CI: 11,800–50,800)
    # after adjusting for case underascertainment and underreporting.
    confidence_interval=(11_800, 50_800),  # 95% Bootstrap Confidence Interval
    coverage_probability=0.95,
    country="United States",
    # I picked 2019 rates, as estimated rates in 2020 were 30% lower, even
    # while deaths stayed the same, pointing to undercounting due to
    # COVID-19: https://www.cdc.gov/hepatitis/statistics/2020surveillance/introduction/national-profile.htm#:~:text=hepatitis%20B%20transmission.-,Data%20from%20death%20certificates%20filed%20in%20the%20vital%20records%20offices%20of,same%20as%20the%20rate%20during%202019%20(0.42%20deaths%20per%20100%2C000%20population).,-Hepatitis%20C
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


def estimate_prevalences() -> list[Prevalence]:
    chronic_2020 = estimated_chronic_us_2020.to_rate(us_population(year=2020))
    # Hep B chronic cases should be approximately constant, so we
    # can use chronic 2020 estimates for 2021.
    chronic_2021 = dataclasses.replace(
        chronic_2020, date_source=Variable(date="2021")
    )
    dk_2015 = dataclasses.replace(
        denmark_estimated_chronic_2007, date_source=Variable(date="2015")
    )
    dk_2016 = dataclasses.replace(
        denmark_estimated_chronic_2007, date_source=Variable(date="2016")
    )
    dk_2017 = dataclasses.replace(
        denmark_estimated_chronic_2007, date_source=Variable(date="2017")
    )
    dk_2018 = dataclasses.replace(
        denmark_estimated_chronic_2007, date_source=Variable(date="2018")
    )

    return [
        denmark_estimated_chronic_2007,
        dk_2015,
        dk_2016,
        dk_2017,
        dk_2018,
        chronic_2020,
        chronic_2021,
    ]


def estimate_incidences() -> list[IncidenceRate]:
    # Before the COVID-19 pandemic, acute HBV incidence has stayed
    # approximately constant: (https://www.cdc.gov/hepatitis/statistics/
    # 2019surveillance/Figure2.1.htm). We can, therefore, estimate 2020
    # and 2021 incidence from 2019 incidence. We don't expect Hep B to be
    # affected by COVID-19 itself as strongly, given that it is transmitted
    # through sexual contact and injection drug use, which are less affected
    # by social distancing.
    acute_2019 = estimated_acute_2019.to_rate(us_population_2019)

    acute_2020 = dataclasses.replace(
        acute_2019, date_source=Variable(date="2020")
    )
    acute_2021 = dataclasses.replace(
        acute_2019, date_source=Variable(date="2021")
    )
    return []
