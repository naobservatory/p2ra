import dataclasses

from pathogen_properties import *

background = """Herpes Simplex Virus 1 is a very common herpesvirus that causes oral herpes (CDC, https://www.
hopkinsmedicine.org/health/conditions-and-diseases/herpes-hsv1-and-hsv2), manifesting as cold sores or fever
blisters (https://www.hopkinsmedicine.org/health/conditions-and-diseases/herpes-hsv1-and-hsv2). After initial
infection and potential symptoms, most HSV-1 infections persist lifelong."""


pathogen_chars = PathogenChars(
    na_type=NAType.DNA,
    enveloped=Enveloped.ENVELOPED,
    taxid=TaxID(10298),
    selection=SelectionRound.ROUND_1,
)


cdc_2015_2016_nhanes_estimate = Prevalence(
    infections_per_100k=0.478 * 100_000,
    confidence_interval=(0.4281 * 100_000, 0.5277 * 100_000),
    # This CDC estimate is based on the following NHANES data:
    # https://web.archive.org/web/20220707050306/https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/HSV_I.htm#:~:text=YEARS%20%2D%2049%20YEARS-,LBXHE1%20%2D%20Herpes%20Simplex%20Virus%20Type%201,-Variable%20Name%3A
    country="United States",
    start_date="2015",
    end_date="2016",
    active=Active.LATENT,
    source="https://www.cdc.gov/nchs/data/databriefs/db304_table.pdf#page=1",
    methods="https://www.cdc.gov/nchs/products/databriefs/db304.htm#:~:text=Data%20for%20this,p%20%3C%200.05",
)

german_seroprevalence_2008_2011 = Prevalence(
    infections_per_100k=0.787 * 100_000,
    confidence_interval=(0.772 * 100_000, 0.801 * 100_000),  # 95% CI
    # "Weighted seroprevalence of HSV1:
    # The overall seroprevalence of HSV1 in the DEGS was 78.7% (95%CI 77.2–80.
    # 1).""
    number_of_participants=6627,
    # "In total, 6627 DEGS participants were tested for HSV1, and 5013 were
    # also tested for HSV2.""
    country="Germany",
    start_date="2008",
    end_date="2011",
    active=Active.LATENT,
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5500947/#:~:text=Weighted%20seroprevalence%20of,CI%2077.2%E2%80%9380.1).",
)


def estimate_prevalences() -> list[Prevalence]:
    # HSV_1 prevalence should be close to constant, so extrapolate from
    # 2015-2016 to 2020 and 2021.
    us_2020 = dataclasses.replace(
        cdc_2015_2016_nhanes_estimate, date_source=Variable(date="2020")
    )
    us_2021 = dataclasses.replace(
        cdc_2015_2016_nhanes_estimate, date_source=Variable(date="2021")
    )

    # We assume that the demographics of Denmark and Germany are similar
    # enough to extrapolate German seroprevalence data to Denmark.
    # We furthermore assume that HSV-1 prevalence remains constant over time,
    # extrapolating measurements to 2015-2018.

    dk_2015 = dataclasses.replace(
        german_seroprevalence_2008_2011,
        date_source=Variable(date="2015"),
        location_source=Variable(country="Denmark"),
    )
    dk_2016 = dataclasses.replace(
        german_seroprevalence_2008_2011,
        date_source=Variable(date="2016"),
        location_source=Variable(country="Denmark"),
    )
    dk_2017 = dataclasses.replace(
        german_seroprevalence_2008_2011,
        date_source=Variable(date="2017"),
        location_source=Variable(country="Denmark"),
    )
    dk_2018 = dataclasses.replace(
        german_seroprevalence_2008_2011,
        date_source=Variable(date="2018"),
        location_source=Variable(country="Denmark"),
    )

    return [
        us_2020,
        us_2021,
        dk_2015,
        dk_2016,
        dk_2017,
        dk_2018,
    ]


def estimate_incidences() -> list[IncidenceRate]:
    return []
