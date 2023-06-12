from pathogen_properties import *

background = """Hepatitis C is a chronic liver condition, caused by the hepatitis C virus. It is most often transmitted through IV drug use, but also through birth and sexual intercourse."""

pathogen_chars = PathogenChars(
    na_type=NAType.RNA,
    enveloped=Enveloped.ENVELOPED,
    taxid=TaxID(11103),
)

us_population_2019 = Population(
    people=328.2 * 1e6,
    country="United States",
    date="2019",
    source="https://data.census.gov/table?q=us+population+2019&tid=ACSDP1Y2019.DP05",
)


# Data below is from the Ohio Department of Health. They give case rates,
# which are "shown per 100,000 persons and were calculated using census
# estimates for that year, except 2021 is using 2020 census."


reported_acute_ohio_2020 = IncidenceRate(
    annual_infections_per_100k=2.0,
    # "Most commonly, acute hepatitis C virus (HCV) infection is defined as
    # the 6-month time period following acquisition of hepatitis C virus.
    # (https://www.hepatitisc.uw.edu/go/screening-diagnosis/acute-diagnosis/core-concept/all#:~:text=Most%20commonly%2C%20acute%20hepatitis%20C%20virus%20(HCV)%20infection%20is%20defined%20as%20the%206%2Dmonth%20time%20period%20following%20acquisition%20of%20hepatitis%20C%20virus.)"
    date="2020",
    country="United States",
    state="Ohio",
    source="https://odh.ohio.gov/wps/wcm/connect/gov/ec0dec22-1eea-4d17-a86a-ac4bc35be4d3/HCV+5+Year+Report+2021.pdf?MOD=AJPERES&CONVERT_TO=url&CACHEID=ROOTWORKSPACE.Z18_M1HGGIK0N0JO00QO9DDDDM3000-ec0dec22-1eea-4d17-a86a-ac4bc35be4d3-oqU9kQ8"
    # Page 1 of 8
)

reported_acute_ohio_2021 = IncidenceRate(
    annual_infections_per_100k=1.3,
    date="2021",
    country="United States",
    state="Ohio",
    source="https://odh.ohio.gov/wps/wcm/connect/gov/ec0dec22-1eea-4d17-a86a-ac4bc35be4d3/HCV+5+Year+Report+2021.pdf?MOD=AJPERES&CONVERT_TO=url&CACHEID=ROOTWORKSPACE.Z18_M1HGGIK0N0JO00QO9DDDDM3000-ec0dec22-1eea-4d17-a86a-ac4bc35be4d3-oqU9kQ8",
)


estimated_chronic_us_2013_2016 = Prevalence(
    # [...] we analyzed 2013-2016 data from the National Health
    # and Nutrition Examination Survey (NHANES) to estimate the prevalence of
    # HCV in the noninstitutionalized civilian population and used a
    # combination of literature reviews and population size estimation
    # approaches to estimate the HCV prevalence and population sizes for four
    # additional populations: incarcerated people, unsheltered homeless
    # people, active-duty military personnel, and nursing home residents.
    infections_per_100k=0.01 * 100_000,  # among all US adults
    # We estimated that [...] 1.0% (95% CI, 0.8-1.1%) of all adults,
    # approximately 2.4 (2.0-2.8) million persons, were HCV RNA–positive
    # (indicating current infection).
    confidence_interval=(0.008 * 100_000, 0.011 * 100_000),  # 95% CI
    coverage_probability=0.95,
    date="2016",
    # To calculate the number of noninstitutionalized civilians in the United
    # States with HCV antibody and HCV RNA during 2013-2016, prevalence
    # estimates were then multiplied by the estimated total
    # noninstitutionalized civilian adult U.S. population as of December 31,
    # 2016 from the 2012-2016 ACS.
    country="United States",
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6719781/#:~:text=(1%2C983%2C900%2D2%2C807%2C800)-,1.0%25,(0.8%25%2D1.1%25),-Open%20in%20a",
    active=Active.LATENT,
)

cdc_estimated_acute_2019 = IncidenceAbsolute(
    annual_infections=57_500,
    # During 2019, a total of 4,136 acute hepatitis C cases were reported to
    # CDC, corresponding to 57,500 estimated infections (95% CI: 45,500–196,
    # 000).
    # We do not use 2020 estimated incidence, as the underascertainment factor
    # with which the CDC gets from reported to estimated cases wasn't changed
    # in 2020, even though the pandemic should have influenced case reporting
    # significantly.
    confidence_interval=(45_500, 196_000),  # 95% CI
    coverage_probability=0.95,
    date="2019",
    country="United States",
    source="https://www.cdc.gov/hepatitis/statistics/2019surveillance/Introduction.htm#Technical:~:text=During%202019%2C%20a%20total%20of%204%2C136%20acute%20hepatitis%20C%20cases%20were%20reported%20to%20CDC%2C%20corresponding%20to%2057%2C500%20estimated%20infections%20(95%25%20CI%3A%2045%2C500%E2%80%93196%2C000))",
)

acute_underreporting_factor = Scalar(
    scalar=13.9,
    confidence_interval=(11.0, 47.4),  # 95% CI
    coverage_probability=0.95,
    # each reported case of acute hepatitis C represents 13.9 estimated
    # infections (95% bootstrap CI: 11.0–47.4).
    source="https://www.cdc.gov/hepatitis/statistics/2018surveillance/pdfs/2018HepSurveillanceRpt.pdf?#page=8",
)


ohio_counties_case_rates = {
    # source: https://odh.ohio.gov/wps/wcm/connect/gov/ec0dec22-1eea-4d17-a86a-ac4bc35be4d3/HCV+5+Year+Report+2021.pdf?MOD=AJPERES&CONVERT_TO=url&CACHEID=ROOTWORKSPACE.Z18_K9I401S01H7F40QBNJU3SO1F56-ec0dec22-1eea-4d17-a86a-ac4bc35be4d3-oqU9kQ8
    # We ended up not using total cases, which could be used to arrive at
    # chronic incidence by subtracting acute cases [total - acute]. This is
    # because because we i) do not have an underreporting factor for chronic
    # incidence, and ii) prevalence is a better measure for the chronic
    # version of Hep C.
    "Franklin": {
        "2020": {"acute": 0.8, "total": 74.3},
        "2021": {"acute": 1.5, "total": 86.8},
    },
    "Greene": {
        "2020": {"acute": 0.0, "total": 74.1},
        "2021": {"acute": 0.6, "total": 49.4},
    },
    "Lawrence": {
        "2020": {"acute": 15.2, "total": 367.2},
        "2021": {"acute": 1.7, "total": 379.1},
    },
    "Licking": {
        "2020": {"acute": 3.4, "total": 59.5},
        "2021": {"acute": 2.2, "total": 65.1},
    },
    "Lucas": {
        "2020": {"acute": 0.7, "total": 134.7},
        "2021": {"acute": 0.2, "total": 122.1},
    },
    "Montgomery": {
        "2020": {"acute": 0.6, "total": 105.9},
        "2021": {"acute": 0.4, "total": 98.0},
    },
    "Sandusky": {
        "2020": {"acute": 1.7, "total": 104.5},
        "2021": {"acute": 0.0, "total": 102.8},
    },
    "Summit": {
        "2020": {"acute": 1.1, "total": 99.7},
        "2021": {"acute": 0.7, "total": 97.2},
    },
    "Trumbull": {
        "2020": {"acute": 6.1, "total": 120.9},
        "2021": {"acute": 1.0, "total": 107.2},
    },
}

OHIO_COUNTY_ESTIMATES_SOURCE = "https://odh.ohio.gov/wps/wcm/connect/gov/ec0dec22-1eea-4d17-a86a-ac4bc35be4d3/HCV+5+Year+Report+2021.pdf?MOD=AJPERES&CONVERT_TO=url&CACHEID=ROOTWORKSPACE.Z18_K9I401S01H7F40QBNJU3SO1F56-ec0dec22-1eea-4d17-a86a-ac4bc35be4d3-oqU9kQ8"


def estimate_incidences():
    # Hep C acute cases should be approximately constant, so we
    # can use 2019 acute estimates for 2020 and 2021.

    acute_2019 = cdc_estimated_acute_2019.to_rate(us_population_2019)
    acute_2020 = dataclasses.replace(
        acute_2019, date_source=Variable(date="2020")
    )
    acute_2021 = dataclasses.replace(
        acute_2019, date_source=Variable(date="2021")
    )

    estimates = [
        acute_2019,
        acute_2020,
        acute_2021,
        reported_acute_ohio_2020 * (acute_underreporting_factor),
        reported_acute_ohio_2021 * (acute_underreporting_factor),
    ]
    for county in ohio_counties_case_rates:
        for year in ohio_counties_case_rates[county]:
            estimates.append(
                IncidenceRate(
                    annual_infections_per_100k=ohio_counties_case_rates[
                        county
                    ][year]["acute"],
                    date=year,
                    country="United States",
                    state="Ohio",
                    county=county,
                    source=OHIO_COUNTY_ESTIMATES_SOURCE,
                )
                * acute_underreporting_factor
            )
    return estimates


def estimate_prevalences() -> list[Prevalence]:
    # This estimate being from 2016, I (Simon) do not want to transfer it to
    # 2019-202 without first doing further research.
    estimates = [estimated_chronic_us_2013_2016]
    return estimates
