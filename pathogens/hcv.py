from pathogen_properties import *

background = """Hepatitis C is a chronic liver condition, caused by the hepatitis C virus. It is most often transmitted through IV drug use, but also through birth and sexual intercourse."""


# Get state_level numbers for Ohio!

pathogen_chars = PathogenChars(
    na_type=NAType.RNA,
    enveloped=Enveloped.ENVELOPED,
    taxid=TaxID(11103),
)

# Data below is for Ohio only, targeting the Spurbeck paper (https://www.frontiersin.org/articles/10.3389/fpubh.2023.1145275/full).
# Data for all other states can be found here: https://www.cdc.gov/hepatitis/statistics/2020surveillance/hepatitis-c/table-3.1.htm

reported_acute_ohio_2020 = IncidenceRate(
    # The Ohio Department of Health:
    # "Rates are shown per 100,000 persons and were calculated using census estimates for that year, except 2021 is using 2020 census.""
    annual_infections_per_100k=2,
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
    date=2020,
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
)
