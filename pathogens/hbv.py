from pathogen_properties import *

hbv_background = """Hepatitis B is a liver infection caused by the hepatitis B virus. It is transmitted through birth and contact with infected blood or bodily fluids. Hepatitis B infection increases the risk for hepatocellular carcinoma"""

# TODO: Add case surveillance data: https://wonder.cdc.gov/nndss/static/2019/annual/2019-table2h.html


pathogen_chars = PathogenChars(
    na_type=NAType.DNA,
    enveloped=Enveloped.ENVELOPED,
    taxid=TaxID(10407),
)

#To add: in-house NHANES estimate

Prevalence(
    infections_per_100k=
)
