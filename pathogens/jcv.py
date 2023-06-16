import dataclasses

from pathogen_properties import *

background = """JC virus is a common virus which has minimal impact in
immunocompetent humans"""

pathogen_chars = PathogenChars(
    na_type=NAType.DNA,
    enveloped=Enveloped.NON_ENVELOPED,
    taxid=TaxID(10632),
)


def estimate_prevalences() -> list[Prevalence]:
    stub_estimate = Prevalence(
        # Completely made up, checking in stub for model testing.
        infections_per_100k=50_000,
        active=Active.LATENT,
        country="United States",
    )
    return [
        # This pathogen should be close to constant, so extrapolate to 2020 and
        # 2021.
        dataclasses.replace(stub_estimate, date_source=Variable(date="2020")),
        dataclasses.replace(stub_estimate, date_source=Variable(date="2021")),
    ]


def estimate_incidences():
    return []
