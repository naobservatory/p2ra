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
    # Among 15–59-year-olds, 2013–2016 prevalence of any HPV infection was 40%
    # overall, 41.8% among males, and 38.4% among females.
    infections_per_100k=0.400 * 100_000,
    # We assume that measurements for 18 to 59 year olds roughly correspond to
    # the true all-age population prevalence.
    confidence_interval=(0.392 * 100_000, 0.409 * 100_000),
    coverage_probability=0.5,  # confidence_interval
    number_of_participants=8005,
    country="United States",
    start_date="2013",
    end_date="2016",
    active=Active.LATENT,
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10037549/#:~:text=Any%20HPV-,40.0%20(39.2%2C%2040.9),-41.8%20(40.6%2C%2042.9",
)

nhanes_based_2018_18_59_yo_incidence = IncidenceRate(
    annual_infections_per_100k=1222 * 10,  # incidence per 10k -> per 100k.
    # We assume that measurements for 18 to 59 year olds roughly correspond to
    # the true all-age population prevalence.
    confidence_interval=(969 * 10, 1436 * 10),
    coverage_probability=0.5,  # uncertainty interval
    country="United States",
    date="2018",
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10037549/#:~:text=Any%20HPV-,1222%20(969%2C%201436),-1223%20(983%2C%201332",
)


def estimate_prevalences():
    # HPV prevalence is probably decreasing from 2013-2016 to the present, as
    # vaccination is increasing and they're covering more subtypes.
    # Extrapolating to 2020 and 2021 without adjusting for this might not be
    # very good.
    # TODO(#157): look into whether we can get estimates for 2020 and 2021.
    return []


def estimate_incidences():
    # HPV prevalence should be close to constant, so extrapolate from
    # 2018 to 2020 and 2021.
    return [
        dataclasses.replace(
            nhanes_based_2018_18_59_yo_incidence,
            date_source=Variable(date="2020"),
        ),
        dataclasses.replace(
            nhanes_based_2018_18_59_yo_incidence,
            date_source=Variable(date="2021"),
        ),
    ]
