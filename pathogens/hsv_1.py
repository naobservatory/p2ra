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

tear_and_saliva_prevalence = Prevalence(
    infections_per_100k=0.98 * 100_000,
    number_of_participants=50,
    country="United States",
    state="Louisiana",
    date="2005",
    active=Active.LATENT,
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1200985/#:~:text=Frequency%20of%20HSV%2D1%20DNA%20shedding%20over,swab%20data%20(subjects%2030%2C%2038%20not%20available).",
)


def estimate_prevalences() -> list[Prevalence]:
    return [
        # HSV_1 prevalence should be close to constant, so extrapolate from
        # 2015-2016 to 2020 and 2021.
        dataclasses.replace(
            cdc_2015_2016_nhanes_estimate, date_source=Variable(date="2020")
        ),
        dataclasses.replace(
            cdc_2015_2016_nhanes_estimate, date_source=Variable(date="2021")
        ),
    ]


def estimate_incidences() -> list[IncidenceRate]:
    return []
