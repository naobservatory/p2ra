from pathogen_properties import *

background = """Dengue virus, a mosquito-borne viral infection primarily transmitted by the Aedes aegypti and Aedes albopictus mosquitoes, present in tropical and subtropical regions. The virus presents with flu-like symptoms that can progress to severe dengue, characterized by hemorrhagic fever, shock, and even death."""


pathogen_chars = PathogenChars(
    na_type=NAType.RNA,
    enveloped=Enveloped.ENVELOPED,
    taxid=12814,
)

california_reported_cases = PrevalenceAbsolute(
    infections=86,
    country="United States",
    state="California",
    date="2020",
    source="https://www.cdph.ca.gov/Programs/CID/DCDC/CDPH%20Document%20Library/TravelAssociatedCasesofDengueVirusinCA.pdf",
    tag="California-2020",
)

CA_population = Population(
    people=39000000,
    tag="California-2020",
)
# disease duration in days
disease_duration = Scalar(
    scalar=7,
    source="https://www.cdc.gov/dengue/symptoms/index.html",
)

per100k = Scalar(
    scalar=100000,
)

reciprocal_days_in_a_year = Scalar(
    scalar=1 / 365,
)


def estimate_prevalences():
    return [
        california_reported_cases.to_rate(CA_population)
        .scale(disease_duration)
        .scale(per100k)
        .scale(reciprocal_days_in_a_year)
    ]