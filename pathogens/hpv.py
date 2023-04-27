from pathogen_properties import *

background = """Human papilloma virus is a species of virus that has over 100 different subtypes that have varying clinical significance. Together HPV-infections are among the most common STIs globally, though prevalence is in decline following vaccination drives. Vaccines mostly target high-risk HPV subtypes that are associated with cancer."""

#TODO: Discuss shedding duration -> should we use mean or median?

pathogen_chars = PathogenChars(
    na_type=NAType.RNA,
    enveloped=Enveloped.NON_ENVELOPED,
    taxid=10566,
)

median_shedding_duration_brazil_1993_1997 = SheddingDuration(
    days=8.2 * 30,
    confidence_interval=(8.0 * 30,9.2 * 30),
    country="Brazil",
    city="Sao Paulo",
    start_date="1993",
    end_date="1997",
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7889327/#:~:text=For%20any%20HPV,P%20%3D%20.009).",
    )



cdc_2013_2016_nhanes_estimate = Prevalence(
    infections_per_100k=0.400 * 100_000,
    confidence_interval=(0.392 * 100_000, 0.409 * 100_000),
    country="United States",
    start_date="2013",
    end_date="2016",
    source
    methods=
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10037549/#:~:text=Any%20HPV-,40.0%20(39.2%2C%2040.9),-41.8%20(40.6%2C%2042.9"







def estimate_prevalences():
    return [
        us_incidence_absolute_2018.to_rate(us_population_2018).to_prevalence(
            hva_shedding_duration
        ),
        king_county_confirmed_cases_rate_2017.to_prevalence(
            hva_shedding_duration
        ).scale(incidence_underreporting_scalar),
        king_county_confirmed_cases_rate_2018.to_prevalence(
            hva_shedding_duration
        ).scale(incidence_underreporting_scalar),
    ]
