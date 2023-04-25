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


cdc_2015_2016_nhanes_seroprevalence = Prevalence(
    infections_per_100k=0.496 * 100_000,
    number_of_participants=3710,
    start_date="2015",
    end_date="2016",
    country="United States",
    source="https://web.archive.org/web/20220707050306/https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/HSV_I.htm",
)
cdc_2015_2016_nhanes_estimate = Prevalence(
    infections_per_100k=0.478 * 100_000,
    confidence_interval=(0.4281 * 100_000, 0.5277 * 100_000),
    country="United States",
    start_date="2015",
    end_date="2016",
    source="https://www.cdc.gov/nchs/data/databriefs/db304_table.pdf#page=1",
    methods="https://www.cdc.gov/nchs/products/databriefs/db304.htm#:~:text=Data%20for%20this,p%20%3C%200.05",
)

tear_and_saliva_prevalence = Prevalence(
    infections_per_100k=0.98 * 100_000,
    number_of_participants=50,
    country="United States",
    state="Louisiana",
    date="2005",
    source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1200985/#:~:text=Frequency%20of%20HSV%2D1%20DNA%20shedding%20over,swab%20data%20(subjects%2030%2C%2038%20not%20available).",
)


def estimate_prevalences():
    return [
        cdc_2015_2016_nhanes_seroprevalence,
        cdc_2015_2016_nhanes_estimate,
        tear_and_saliva_prevalence,
    ]
