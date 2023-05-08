from pathogen_properties import *

background = """Varicella-zoster virus (VZV) causes both varicella (chickenpox) and herpes zoster (shingles). After a primary infection with VZV,
the virus can remain dormant in the sensory nerve cells of the body,
reactivating later in life and causing herpes zoster, or shingles."""


pathogen_chars = PathogenChars(
    na_type=NAType.DNA,
    enveloped=Enveloped.ENVELOPED,
    taxid=TaxID(10335),
)

#TODO Either delete or turn the data below into national prevalence
# measurements, identifying an appropriate adjustment factor.

us_2022_provisional_varicella_cases = IncidenceAbsolute(
    date="2022",
    #Files containing weekly cases can be downloaded with the bash script
    # download_vzv_cases.sh, found under ../prevalence-data/ . Given the low
    # case weekly case counts, we decided to not use that data.
    country="US",
    annual_infections=3890,
    source="https://wonder.cdc.gov/nndss/static/2022/52/2022-52-table1kk.html#:~:text=119-,3%2C890,-3%2C665"
    #These are provisional cases from the National Notifiable Diseases Surveillance System (NNDSS) by the CDC, taken from the cumulate count in of the final weekly report of 2022.
)

us_2021_provisional_varicella_cases = IncidenceAbsolute(
    date="2021",
    country="US",
    annual_infections=3,665,
    source="https://wonder.cdc.gov/nndss/static/2022/52/2022-52-table1kk.html#:~:text=3%2C890-,3%2C665,-New%20England"
    #Same disclaimer as above.
)


us_2019_reported_varicella_cases = IncidenceAbsolute(
    date="2019",
    country="US",
    annual_infections=8297,
    source="https://wonder.cdc.gov/nndss/static/2019/annual/2019-table2r.html#:~:text=3-,8%2C297,-6"
    #These are yearly reported cases from the National Notifiable Diseases Surveillance System (NNDSS) by the CDC. Case counts are likely more complete and thus higher than the provisional counts.

us_2011_shingles_incidence_insured = IncidenceAbsolute(
    # QUESTION: We could scale this number to the US population by accounting
    # for the participant's age distribution, viewable here: https://www.ncbi.
    # nlm.nih.gov/pmc/articles/PMC4636742/bin/12879_2015_1262_MOESM1_ESM.doc
    date="2011",
    country="US",
    annual_infections=113_861,
    number_of_participants=27_616_373,
    # Note that this sample is skewed toward the commercially insured.
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4636742/#:~:text=N%E2%80%89%3D-,113%2C%20861,-N%E2%80%89%3D%E2%80%8927%2C616%2C373"
)


# Notes on vaccination

"""VZV vaccination among seniors looks very effective:

Our analysis included 3·36 million person-years of data, corresponding to an average of 310 001 patients aged 60–89 years who were registered at an RCGP practice each year. By Aug 31, 2016, uptake of the vaccine varied between 58% for the recently targeted cohorts and 72% for the first routine cohort. Across the first 3 years of vaccination for the three routine cohorts, incidence of herpes zoster fell by 35% (incidence rate ratio 0·65 [95% 0·60–0·72]) and of postherpetic neuralgia fell by 50% (0·50 [0·38–0·67]). The equivalent reduction for the four catch-up cohorts was 33% for herpes zoster (incidence rate ratio 0·67 [0·61–0·74]) and 38% for postherpetic neuralgia (0·62 [0·50–0·79]). These reductions are consistent with a vaccine effectiveness of about 62% against herpes zoster and 70–88% against postherpetic neuralgia.
"""

"""
Shingles vaccination rate has increased from 6.7 to 34.5 in between 2008 to 2018 (https://www.cdc.gov/nchs/data/databriefs/db370-tables-508.pdf#page=1)
"""


#TODO: Do more research here
"""
Nevertheless, Herpes Zoster incidence has icnreased continuously, adjusting for age. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4928389/#:~:text=Of%20the%208017,varicella%20vaccination%20program."


shingles_duration = SheddingDuration(
    days=28,
    source="https://www.nia.nih.gov/health/shingles#:~:text=Most%20cases%20of%20shingles%20last%20three%20to%20five%20weeks.",
)

chickenpox_duration = SheddingDuration(
    days=12,
    source="https://my.clevelandclinic.org/health/diseases/4017-chickenpox#:~:text=Chickenpox%20usually%20goes%20away%20after%2010%20to%2014%20days.",
)
