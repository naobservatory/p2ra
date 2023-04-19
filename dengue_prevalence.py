#!/usr/bin/env python3

from pathogen_properties import *
import math

background = """Dengue virus, a mosquito-borne viral infection primarily transmitted by the Aedes aegypti and Aedes albopictus mosquitoes, present in tropical and subtropical regions. The virus presents with flu-like symptoms that can progress to severe dengue, characterized by hemorrhagic fever, shock, and even death."""


pathogen_chars = PathogenChars(
    na_type=NAType.RNA.value,
    enveloped=Enveloped.ENVELOPED.value,
    taxid=12814,
)


prevalence_vars = {
    "california_reported_cases_estimate": PrevalenceVariable(
        variable_type=VariableType.EXTERNAL_ESTIMATE.value,
        percentage=0.000005,
        number_of_participants=40000000,
        country="United States",
        state = "California",
        source="https://www.cdph.ca.gov/Programs/CID/DCDC/CDPH%20Document%20Library/TravelAssociatedCasesofDengueVirusinCA.pdf",
    ),
}
