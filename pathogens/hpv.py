from pathogen_properties import *
from populations import us_population


background = """Human papilloma virus is DNA virus with over 100 different 
subtypes that have varying clinical significance. Together HPV-infections are 
among the most common STIs globally, though prevalence is in decline following 
vaccination drives. Vaccines mostly target high-risk HPV subtypes that are 
associated with cancer."""

# TODO: Incorporate NHANES 2017-2018 adult female data:
# https://wwwn.cdc.gov/Nchs/Nhanes/limited_access/HPVW_J_R.htm
# TODO: Incorprate NHANES 2015-2015 male and female genital and oral HPV data:
# https://wwwn.cdc.gov/nchs/nhanes/search/datapage.aspx?Component=Laboratory&CycleBeginYear=2015#:~:text=September%202017-,Human%20Papillomavirus%20(HPV)%20%2D%20Oral%20Rinse,November%202018,-Insulin

pathogen_chars = PathogenChars(
    na_type=NAType.DNA,
    enveloped=Enveloped.NON_ENVELOPED,
    taxid=10566,
)


cdc_2013_2016_nhanes_prevalence_estimate = Prevalence(
    infections_per_100k=0.400 * 100_000,
    confidence_interval=(0.392 * 100_000, 0.409 * 100_000),
    coverage_probability=0.5,  # confidence_interval
    number_of_participants=8005,
    country="United States",
    start_date="2013",
    end_date="2016",
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10037549/#:~:text=Any%20HPV-,40.0%20(39.2%2C%2040.9),-41.8%20(40.6%2C%2042.9",
)


def estimate_prevalences():
    return [cdc_2013_2016_nhanes_prevalence_estimate]


def estimate_incidences():
    return []
