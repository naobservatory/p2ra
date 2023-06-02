from pathogen_properties import *
from populations import us_population

background = """Dengue virus, a mosquito-borne viral infection primarily
transmitted by the Aedes aegypti and Aedes albopictus mosquitoes, present 
in tropical and subtropical regions. The virus presents with flu-like symptoms
that can progress to severe dengue, characterized by
hemorrhagic fever, shock, and even death."""


pathogen_chars = PathogenChars(
    na_type=NAType.RNA,
    enveloped=Enveloped.ENVELOPED,
    taxid=TaxID(12637),
)

california_reported_cases = IncidenceAbsolute(
    annual_infections=86,
    country="United States",
    state="California",
    date="2020",
    source="https://www.cdph.ca.gov/Programs/CID/DCDC/CDPH%20Document%20Library/TravelAssociatedCasesofDengueVirusinCA.pdf",
)


def estimate_incidences() -> list[IncidenceRate]:
    return [
        california_reported_cases.to_rate(
            us_population(state="California", year=2020)
        )
    ]


def estimate_prevalences() -> list[Prevalence]:
    return []
