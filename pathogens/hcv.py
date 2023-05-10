from pathogen_properties import *

background = """Hepatitis C is a chronic liver condition, caused by the hepatitis C virus. It is most often transmitted through IV drug use, but also through birth and sexual intercourse."""


# Get state_level numbers for Ohio!

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


reported_chronic_ohio_2020 = IncidenceRate(
    # Chronic Hepatitis C Virus (HCV) infection is commonly defined by a
    # positive test for HCV RNA 6 months after initial diagnosis of HCV.
    annual_infections_per_100k=110.6,
    date="2020",
    country="United States",
    state="Ohio",
    source="https://odh.ohio.gov/wps/wcm/connect/gov/ec0dec22-1eea-4d17-a86a-ac4bc35be4d3/HCV+5+Year+Report+2021.pdf?MOD=AJPERES&CONVERT_TO=url&CACHEID=ROOTWORKSPACE.Z18_M1HGGIK0N0JO00QO9DDDDM3000-ec0dec22-1eea-4d17-a86a-ac4bc35be4d3-oqU9kQ8",
)


reported_chronic_ohio_2021 = IncidenceRate(
    annual_infections_per_100k=105.4,
    date="2021",
    country="United States",
    state="Ohio",
    source="https://odh.ohio.gov/wps/wcm/connect/gov/ec0dec22-1eea-4d17-a86a-ac4bc35be4d3/HCV+5+Year+Report+2021.pdf?MOD=AJPERES&CONVERT_TO=url&CACHEID=ROOTWORKSPACE.Z18_M1HGGIK0N0JO00QO9DDDDM3000-ec0dec22-1eea-4d17-a86a-ac4bc35be4d3-oqU9kQ8",
)


cdc_estimated_acute_2020 = IncidenceAbsolute(
    # corresponding to 66,700 estimated infections (95% CI: 52,700–227,400) #after adjusting for case underascertainment and underreporting
    annual_infections=66_700,
    confidence_interval=(52_700, 227_400),  # 95% CI
    coverage_probability=0.95,
    date="2020",
    country="United States",
    tag="us-2020",
)
estimated_chronic_us_2013_2016 = PrevalenceAbsolute(
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
    # (indicating current infection)
    confidence_interval=(2_000_000, 2_800_000),  # 95% CI
    coverage_probability=0.95,
    start_date="2013",
    end_date="2016",
    country="United States",
    tag="us-2013-2016",
    source="10.1002/hep.30297",
    active=Active.LATENT
)

SheddingDuration(
    days=3*30.4,  # 3 months
    source="https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0122232#:~:text=Median%20HCV%20RNA%20levels%20initially%20diverged%20at%20three%20months%20following%20infection"
)


ohio_counties_case_counts = {
    #source: https://odh.ohio.gov/wps/wcm/connect/gov/ec0dec22-1eea-4d17-a86a-ac4bc35be4d3/HCV+5+Year+Report+2021.pdf?MOD=AJPERES&CONVERT_TO=url&CACHEID=ROOTWORKSPACE.Z18_K9I401S01H7F40QBNJU3SO1F56-ec0dec22-1eea-4d17-a86a-ac4bc35be4d3-oqU9kQ8
    #county: [acute2020, chronic2020, acute2021, chronic2022]
    "Summit": [1.1, 99.7, 0.7, 97.2],
    "Trumbull": [6.1, 120.9, 1.0, 107.2],
    "Lucas": [0.7, 134.7, 0.2, 122.1],
    "Lawrence": [15.2, 367.2, 1.7, 379.1],
    "Sandusky": [1.7, 104.5, 0.0, 102.8],
    "Franklin": [0.8,74.3,1.5,86.8],
    "Licking": [3.4, 59.5, 2.2, 65.1],
    "Greene": [0.0, 74.1, 0.6, 49.4],
    "Montgomery": [0.6, 105.9, 0.4, 98.0]
}


def estimate_prevalences():
    years = [2020, 2020, 2021, 2021]
    ohio_county_estimates = []
    for key in ohio_counties_case_counts:
        counter = -1
        for rate in ohio_counties_case_counts[key]:
            counter += 1
            if counter == 0 or counter == 2:
                case_rate = IncidenceRate(
                    annual_infections_per_100k=rate,
                    date=years[counter],
                    country="United States",
                    state="Ohio",
                    county=key,
                )
                ohio_county_estimates.append(
                case_rate.to_prevalence(shedding_duration)
            )
            elif counter == 1 or counter == 3:
                prevalence = Prevalence(
                    infections_per_100k=rate,
                    date=years[counter],
                    country="United States",
                    state="Ohio",
                    county=key,
                )
                ohio_county_estimates.append(prevalence)
    return [
        estimates,
        reported_acute_ohio_2020,
        reported_acute_ohio_2021,
        reported_chronic_ohio_2020,
        reported_chronic_ohio_2021,
        cdc_estimated_acute_2020,
        estimated_chronic_us_2013_2016]
