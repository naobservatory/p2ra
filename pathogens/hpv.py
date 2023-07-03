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
    selection=SelectionRound.ROUND_1,
)

dk_male_conscripts_seroprevalence_2016 = Prevalence(
    infections_per_100k=0.418 * 100_000,
    confidence_interval=(0.399 * 100_000, 0.438 * 100_000),  # 95% CI
    # "The overall HPV prevalence was [...] 41.8% (95% CI, 39.9–43.8) with
    # PCR. "
    number_of_participants=2436,
    # "This cross-sectional study, the DanMale study, was conducted between
    # February 2006 and January 2007. Male employees and conscripts at
    # military barracks all over Denmark were invited to participate. [...]
    # The study participants were 18 to 65 years of age; 69.0% were younger
    # than 21 years (mean, 23 years) and 14.5% were 30 years or older."
    country="Denmark",
    date="2016",
    active=Active.LATENT,
    source="https://journals.lww.com/stdjournal/Fulltext/2015/08000/Human_Papillomavirus_Infection_Among_2460_Men_in.12.aspx#:~:text=With%20the%20PCR%20test%2C%20the%20HPV%20prevalence%20was%2041.8%25%3B%20730%20(30.0%25)%20men%20had%20HR%20types%20and%20279%20(11.5%25)%20had%20LR%20types.%20In%20183%20PCR%2Dpositive%20men%20(7.5%25)%2C%20the%20HPV%20type%20could%20not%20be%20identified%20and%20they%20were%20classified%20as%20HPVX%20positive%20(Table%202).",
)

nhanes_2013_2016_18_59_yo_prevalence = Prevalence(
    # Among 15–59-year-olds, 2013–2016 prevalence of any HPV infection was 40%
    # overall, 41.8% among males, and 38.4% among females. These prevalences
    # correspond to 77.3 million persons overall, or 40.5 million males and
    # 37.0 million females, with a prevalent HPV infection in 2018.
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
    # HPV prevalence (any subtype) was extrapolated in the study from
    # 2013-2016 to 2018. We extrapolate forward to 2020 and 2021. Though there
    # is a 9-valent vaccine that got rolled out in 2016, this should not have
    # a very large effect on prevalence, as only young individuals are
    # vaccinated.
    us_2020 = dataclasses.replace(
        nhanes_2013_2016_18_59_yo_prevalence,
        date_source=Variable(date="2020"),
    )
    us_2021 = dataclasses.replace(
        nhanes_2013_2016_18_59_yo_prevalence,
        date_source=Variable(date="2021"),
    )

    # Prevalence among male Danish men matches US estimates. Furthermore,
    # NHANES HPV measurements of men (41.8%) and women (38.4%) are similar,
    # we can thus extrapolate male Danish prevalence to the whole population.
    dk_2015 = dataclasses.replace(
        dk_male_conscripts_seroprevalence_2016,
        date_source=Variable(date="2015"),
    )
    dk_2017 = dataclasses.replace(
        dk_male_conscripts_seroprevalence_2016,
        date_source=Variable(date="2017"),
    )
    dk_2018 = dataclasses.replace(
        dk_male_conscripts_seroprevalence_2016,
        date_source=Variable(date="2018"),
    )
    return [
        us_2020,
        us_2021,
        dk_2015,
        dk_male_conscripts_seroprevalence_2016,
        dk_2017,
        dk_2018,
    ]


def estimate_incidences():
    return []
