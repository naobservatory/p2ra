#!/usr/bin/env python3

from pathogen_properties import *
import math

background = """Herpes Simplex Virus 1 is a very common herpesvirus that causes oral herpes (CDC, https://www.
hopkinsmedicine.org/health/conditions-and-diseases/herpes-hsv1-and-hsv2), manifesting as cold sores or fever
blisters (https://www.hopkinsmedicine.org/health/conditions-and-diseases/herpes-hsv1-and-hsv2). After initial
infection and potential symptoms, most HSV-1 infections persist lifelong."""


pathogen_chars = PathogenChars(
    na_type=NAType.DNA.value,
    enveloped=Enveloped.NON_ENVELOPED.value,
    taxid=10298,
)


prevalence_vars = {
    "cdc_2015_2016_nhanes_seroprevalence": PrevalenceVariable(
        variable_type=VariableType.MEASUREMENT.value,
        percentage=0.496,
        number_of_participants=3710,
        country="United States",
        source="https://web.archive.org/web/20220707050306/https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/HSV_I.htm",
    ),
    "cdc_2015_2016_nhanes_estimate": PrevalenceVariable(
        variable_type=VariableType.EXTERNAL_ESTIMATE.value,
        percentage=0.478,
        country="United States",
        confidence_interval=(0.4281, 0.5277),
        source="https://www.cdc.gov/nchs/data/databriefs/db304_table.pdf#page=1",
        methods="https://www.cdc.gov/nchs/products/databriefs/db304.htm#:~:text=Data%20for%20this,p%20%3C%200.05",
    ),
    "tear_and_saliva_prevalence": PrevalenceVariable(
        variable_type=VariableType.MEASUREMENT.value,
        percentage=0.98,
        number_of_participants=50,
        country="United States",
        state="Louisiana",
        source="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1200985/#:~:text=Frequency%20of%20HSV%2D1%20DNA%20shedding%20over,swab%20data%20(subjects%2030%2C%2038%20not%20available).",
    ),
}
