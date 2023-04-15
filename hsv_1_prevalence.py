# HSV_1 prevalence calculation
# Written by Simon

# Importing enums and math
from enums import *
import math

# ============================

background = """Herpes Simplex Virus 1 is a very common herpesvirus that causes oral herpes (CDC, https://www.
hopkinsmedicine.org/health/conditions-and-diseases/herpes-hsv1-and-hsv2), manifesting as cold sores or fever
blisters (https://www.hopkinsmedicine.org/health/conditions-and-diseases/herpes-hsv1-and-hsv2). After initial
infection and potential symptoms, most HSV-1 infections persist lifelong."""


pathogen_chars = {
    "pathogen_name": "Herpes Simplex Virus 1",
    "pathogen_acronym": "HSV-1",
    "na_type": NAType.DNA.value,
    "enveloped_non_enveloped": enveloped.non_enveloped.value,
    "ncbi_taxid": 10298,
}


prevalence_vars = {
    # In a the 2015-2016 NHANES survey, the seroprevalence of HSV-1 was 49.6% in sample of 3710 people:
    "prevalence_var_1": {
        "variable_name": "cdc_2015_2016_nhanes_seroprevalence",
        "variable_type": "measurement",
        "percentage": 0.496,
        "number_of_participants": 3710,
        "source": "https://web.archive.org/web/20220707050306/https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/HSV_I.htm",
    },
    # Based on the seroprevalence data above, the CDC estimated the prevalence of HSV-1:
    "prevalence_var_2": {
        "variable_name": "cdc_2015_2016_nhanes_estimate",
        "variable_type": "external_estimate",
        "percentage": 0.478,
        "confidence_interval": (0.4281, 0.5277),
        "source": "https://www.cdc.gov/nchs/data/databriefs/db304_table.pdf#page=1",
        "methods": "https://www.cdc.gov/nchs/products/databriefs/db304.htm#:~:text=Data%20for%20this,p%20%3C%200.05",
    },
    # In a 2015 study, 50 individuals self-collected tear and saliva samples twice daily which was tested for HSV 1
    # DNA. 23 Individuals reported prior herpetic disease The sample alse skewed toward African Americans (n=37; 78%).
    # Across the 50-day study period, 49 out of 50 subjects shed HSV DNA atleast once.
    "prevalence_var_3": {
        "variable_name": "tear_and_saliva_prevalence",
        "variable_type": "measurement",
        "percentage": 0.98,  # percentage prevalence
        "number_of_participants": 50,  # Number of people tested
        "source": "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1200985/#:~:text=Frequency%20of%20HSV%2D1%20DNA%20shedding%20over,swab%20data%20(subjects%2030%2C%2038%20not%20available).",
    },
}

# Making variables and functions available to other scripts
__all__ = ["pathogen_chars", "prevalence_vars", "background"]
