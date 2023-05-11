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
    annual_infections_per_100k=2,
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


estimated_current_infection_us_2013_2016 = PrevalenceAbsolute(
    # [...] we analyzed 2013-2016 data from the National Health
    # and Nutrition Examination Survey (NHANES) to estimate the prevalence of
    # HCV in the noninstitutionalized civilian population and used a
    # combination of literature reviews and population size estimation
    # approaches to estimate the HCV prevalence and population sizes for four
    # additional populations: incarcerated people, unsheltered homeless
    # people, active-duty military personnel, and nursing home residents.
    infections=2_400_000,  # among all US adults
    # We estimated that [...] 1.0% (95% CI, 0.8-1.1%) of all adults,
    # approximately 2.4 (2.0-2.8) million persons, were HCV RNA–positive
    # (indicating current infection).
    confidence_interval=(2_000_000, 2_800_000),  # 95% CI
    coverage_probability=0.95,
    start_date="2013",
    end_date="2016",
    country="United States",
    tag="us-2013-2016",
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6719781/#:~:text=2.1%20million%20persons%20(95%25%20CI%20%3D%201.8%20to%202.5%20million%20persons)%20with%20current%20HCV%20infection%20in%20the%20U.S.%20noninstitutionalized%20civilian%20population.",
    active=Active.LATENT,
)

shedding_duration = SheddingDuration(
    days=3 * 30.4,  # 3 months
    # At month two, median HCV RNA levels remained comparable between
    # individuals with persistent infection (5.4 log IU/mL; IQR: 3.1, 6.4) and
    # spontaneous clearance (4.8 log/IU/mL; IQR: 0.0, 6.0; P = 0.38). Median
    # HCV RNA levels initially diverged at three months following infection,
    # being 4.8 log/IU/mL (IQR: 3.3, 6.0) in individuals with persistent
    # infection compared to 3.2 log/IU/mL (IQR: 0.0, 6.1) in those with
    # spontaneous clearance (P = 0.03).
    source="https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0122232#:~:text=Median%20HCV%20RNA%20levels%20initially%20diverged%20at%20three%20months%20following%20infection",
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
    # each reported case of acute hepatitis C represents 13.9 estimated
    # infections (95% bootstrap CI: 11.0–47.4).
    source="https://www.cdc.gov/hepatitis/statistics/2018surveillance/pdfs/2018HepSurveillanceRpt.pdf?#page=8",
)


ohio_counties_case_counts = {
    # source: https://odh.ohio.gov/wps/wcm/connect/gov/ec0dec22-1eea-4d17-a86a-ac4bc35be4d3/HCV+5+Year+Report+2021.pdf?MOD=AJPERES&CONVERT_TO=url&CACHEID=ROOTWORKSPACE.Z18_K9I401S01H7F40QBNJU3SO1F56-ec0dec22-1eea-4d17-a86a-ac4bc35be4d3-oqU9kQ8
    # county: [acute2020, total2020, acute2021, total2021]
    "Summit": [1.1, 99.7, 0.7, 97.2],
    "Trumbull": [6.1, 120.9, 1.0, 107.2],
    "Lucas": [0.7, 134.7, 0.2, 122.1],
    "Lawrence": [15.2, 367.2, 1.7, 379.1],
    "Sandusky": [1.7, 104.5, 0.0, 102.8],
    "Franklin": [0.8, 74.3, 1.5, 86.8],
    "Licking": [3.4, 59.5, 2.2, 65.1],
    "Greene": [0.0, 74.1, 0.6, 49.4],
    "Montgomery": [0.6, 105.9, 0.4, 98.0],
}


def estimate_prevalences():
    years = ["2020", "2020", "2021", "2021"]
    source = "https://odh.ohio.gov/wps/wcm/connect/gov/ec0dec22-1eea-4d17-a86a-ac4bc35be4d3/HCV+5+Year+Report+2021.pdf?MOD=AJPERES&CONVERT_TO=url&CACHEID=ROOTWORKSPACE.Z18_K9I401S01H7F40QBNJU3SO1F56-ec0dec22-1eea-4d17-a86a-ac4bc35be4d3-oqU9kQ8"
    ohio_county_estimates = []
    for key in ohio_counties_case_counts:
        counter = -1
        for rate in ohio_counties_case_counts[key]:
            counter += 1
            print(counter)
            if counter == 0 or counter == 2:
                case_rate = IncidenceRate(
                    annual_infections_per_100k=ohio_counties_case_counts[key][
                        counter
                    ],
                    date=years[counter],
                    country="United States",
                    state="Ohio",
                    county=key,
                    source=source,
                )
                ohio_county_estimates.append(
                    case_rate.to_prevalence(shedding_duration).__mul__(
                        acute_underreporting_factor
                    )
                )
            elif counter == 1 or counter == 3:
                prevalence = Prevalence(
                    infections_per_100k=rate,
                    date=years[counter],
                    country="United States",
                    state="Ohio",
                    county=key,
                    active=Active.LATENT,
                    source=source,
                )
                ohio_county_estimates.append(prevalence)
    return [
        # Unpack the list of ohio county estimates
        *ohio_county_estimates,
        reported_acute_ohio_2020.to_prevalence(shedding_duration).__mul__(
            acute_underreporting_factor
        ),
        reported_acute_ohio_2021.to_prevalence(shedding_duration).__mul__(
            acute_underreporting_factor
        ),
        reported_total_ohio_2020,
        reported_total_ohio_2021,
        cdc_estimated_acute_44_us_states_2020.to_prevalence(shedding_duration),
        estimated_current_infection_us_2013_2016.to_rate(
            us_adult_population_2016
        ),
    ]
