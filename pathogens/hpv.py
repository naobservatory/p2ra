from pathogen_properties import *
from populations import us_population


background = """Human papilloma virus is DNA virus with over 100 different 
subtypes that have varying clinical significance. Together HPV-infections are 
among the most common STIs globally, though prevalence is in decline following 
vaccination drives. Vaccines mostly target high-risk HPV subtypes that are 
associated with cancer."""

# TODO: Incorporate NHANES 2017-2018 adult female data:
# https://wwwn.cdc.gov/Nchs/Nhanes/limited_access/HPVW_J_R.htm
# TODO: Adjust estimates for missing <18 and >
# TODO: Incorprate NHANES 2015-2015 male and female genital and oral HPV data:
# https://wwwn.cdc.gov/nchs/nhanes/search/datapage.aspx?Component=Laboratory&CycleBeginYear=2015#:~:text=September%202017-,Human%20Papillomavirus%20(HPV)%20%2D%20Oral%20Rinse,November%202018,-Insulin

pathogen_chars = PathogenChars(
    na_type=NAType.DNA,
    enveloped=Enveloped.NON_ENVELOPED,
    taxid=TaxID(10566),
)


nhanes_2013_2016_18_59_yo_prevalence = Prevalence(
    infections_per_100k=0.400 * 100_000,
    # Among 15–59-year-olds, 2013–2016 prevalence of any HPV infection was 40.
    # 0% overall, 41.8% among males, and 38.4% among females.
    confidence_interval=(0.392 * 100_000, 0.409 * 100_000),
    coverage_probability=0.5,  # confidence_interval
    number_of_participants=8005,
    country="United States",
    start_date="2013",
    end_date="2016",
    active=Active.LATENT,
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10037549/#:~:text=Any%20HPV-,40.0%20(39.2%2C%2040.9),-41.8%20(40.6%2C%2042.9",
)

nhanes_based_2018_18_59_yo_incidence = IncidenceAbsolute(
    annual_infections=23.6 * 1_000_000,
    # Among 15–59-year-olds, incidence of any HPV infection was 1222 per 10,
    # 000 persons, with 23.6 million persons acquiring any HPV infection in 2018.
    confidence_interval=(18.7 * 1_000_000, 27.8 * 1_000_000),
    coverage_probability=0.5,  # uncertainty interval
    country="United States",
    date="2018",
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10037549/#:~:text=Any%20HPV-,23.6%20(18.7%2C%2027.8),-11.8%20(9.5%2C%2012.9",
)


population_18_to_59 = Population(
    people=327_167_439 - 73_352_242 - 73_085_935,  # total - u_18 - over_59
    source="https://data.census.gov/table?q=2018+acs&tid=ACSST1Y2018.S0101",
    date="2018",
    country="United States",
)


def estimate_prevalences():
    return [nhanes_2013_2016_18_59_yo_prevalence]


def estimate_incidences():
    return [nhanes_based_2018_18_59_yo_incidence.to_rate(population_18_to_59)]
