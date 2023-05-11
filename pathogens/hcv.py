from pathogen_properties import *

background = """Hepatitis C is a chronic liver condition, caused by the hepatitis C virus. It is most often transmitted through IV drug use, but also through birth and sexual intercourse."""

pathogen_chars = PathogenChars(
    na_type=NAType.RNA,
    enveloped=Enveloped.ENVELOPED,
    taxid=TaxID(11103),
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


reported_total_ohio_2020 = Prevalence(
    # Total hepatitis equals all hepatitis C cases, "acute", "chronic", and
    # "perinatal"
    infections_per_100k=110.6,
    date="2020",
    country="United States",
    state="Ohio",
    source="https://odh.ohio.gov/wps/wcm/connect/gov/ec0dec22-1eea-4d17-a86a-ac4bc35be4d3/HCV+5+Year+Report+2021.pdf?MOD=AJPERES&CONVERT_TO=url&CACHEID=ROOTWORKSPACE.Z18_M1HGGIK0N0JO00QO9DDDDM3000-ec0dec22-1eea-4d17-a86a-ac4bc35be4d3-oqU9kQ8",
    active=Active.LATENT,
)


reported_total_ohio_2021 = Prevalence(
    infections_per_100k=105.4,
    date="2021",
    country="United States",
    state="Ohio",
    source="https://odh.ohio.gov/wps/wcm/connect/gov/ec0dec22-1eea-4d17-a86a-ac4bc35be4d3/HCV+5+Year+Report+2021.pdf?MOD=AJPERES&CONVERT_TO=url&CACHEID=ROOTWORKSPACE.Z18_M1HGGIK0N0JO00QO9DDDDM3000-ec0dec22-1eea-4d17-a86a-ac4bc35be4d3-oqU9kQ8",
    active=Active.LATENT,
)


cdc_estimated_acute_44_us_states_2020 = IncidenceRate(
    annual_infections_per_100k=1.5,
    # The CDC gives a CI for its 66,700 estimated infections (95% CI:
    # 52_700–227_400), but not its rate (1.5 cases per 100k). I adopted their
    # absolute count CI to also get a CI for the rate.
    confidence_interval=(
        (52_700 / 66_700) * 1.5,
        (227_400 / 66_700) * 1.5,
    ),  # 95% CI
    coverage_probability=0.95,
    date="2020",
    country="United States",
    # Specifically, this estimate is for 44 US states, as six states and the
    # District of Columbia did not report data, see here for a list of
    # non-reporters:
    # "https://www.cdc.gov/hepatitis/statistics/2020surveillance/introduction/technical-notes.htm#:~:text=Chronic%20hepatitis%20B-,Acute%20hepatitis%20C,-Chronic%20hepatitis%20C"
    source="https://www.cdc.gov/hepatitis/statistics/2020surveillance/introduction/national-profile.htm#:~:text=During%202020%2C%20a%20total%20of%204%2C798%20acute%20hepatitis%20C%20cases%20were%20reported%20to%20CDC%20from%2044%20states%2C%20corresponding%20to%2066%2C700%20estimated%20infections%20(95%25%20CI%3A%2052%2C700%E2%80%93227%2C400)",
)


estimated_current_infection_us_2013_2016 = Prevalence(
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

shedding_duration = SheddingDuration(
    # In this study, HCV RNA levels in blood from several patient cohorts were
    # aggregated retroactively, identifying different patterns of HCV RNA
    # levels among individuals who cleared HCV, versus those who developed a
    # persistent infection.
    days=3 * 30.4,  # 3 months
    # "At month two, median HCV RNA levels [in blood] remained comparable between individuals with persistent infection (5.4 log IU/mL; IQR: 3.1, 6.4) and spontaneous clearance (4.8 log/IU/mL; IQR: 0.0, 6.0; P = 0.38). Median HCV RNA levels initially diverged at three months following infection, being 4.8 log/IU/mL (IQR: 3.3, 6.0) in individuals with persistent infection compared to 3.2 log/IU/mL (IQR: 0.0, 6.1) in those with spontaneous clearance (P = 0.03).""
    source="https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0122232#:~:text=Median%20HCV%20RNA%20levels%20initially%20diverged%20at%20three%20months%20following%20infection",
    # High HCV viral loads in blood are associated with higher viral loads in
    # rectal fluid, according to this study: https://academic.oup.com/cid/article/64/3/284/2452663#:~:text=Detection%20of%20HCV%20in%20rectal%20fluid%20as%20a%20function%20of%20HCV%20VL%20in%20blood.
)

us_adult_population_2016 = Population(
    people=249_448_772,
    source="https://data.census.gov/table?q=2016+population+us&tid=ACSDP1Y2016.DP05",
    date="2016",
    country="United States",
    tag="us-2013-2016",
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


def estimate_prevalences():
    source = "https://odh.ohio.gov/wps/wcm/connect/gov/ec0dec22-1eea-4d17-a86a-ac4bc35be4d3/HCV+5+Year+Report+2021.pdf?MOD=AJPERES&CONVERT_TO=url&CACHEID=ROOTWORKSPACE.Z18_K9I401S01H7F40QBNJU3SO1F56-ec0dec22-1eea-4d17-a86a-ac4bc35be4d3-oqU9kQ8"
    ohio_county_estimates = []
    for county in ohio_counties_case_rates:
        for year in ohio_counties_case_rates[county]:
            case_rate = IncidenceRate(
                annual_infections_per_100k=ohio_counties_case_rates[county][
                    year
                ]["acute"],
                date=year,
                country="United States",
                state="Ohio",
                county=county,
                source=source,
            )
            ohio_county_estimates.append(
                case_rate.to_prevalence(shedding_duration)
                * (acute_underreporting_factor)
            )
            prevalence = Prevalence(
                infections_per_100k=ohio_counties_case_rates[county][year][
                    "total"
                ],
                date=year,
                country="United States",
                state="Ohio",
                county=county,
                active=Active.LATENT,
                source=source,
            )
            ohio_county_estimates.append(prevalence)
    return [
        # Unpack the list of ohio county estimates
        *ohio_county_estimates,
        reported_acute_ohio_2020.to_prevalence(shedding_duration)
        * (acute_underreporting_factor),
        reported_acute_ohio_2021.to_prevalence(shedding_duration)
        * (acute_underreporting_factor),
        reported_total_ohio_2020,
        reported_total_ohio_2021,
        cdc_estimated_acute_44_us_states_2020.to_prevalence(shedding_duration),
        estimated_current_infection_us_2013_2016,
    ]
