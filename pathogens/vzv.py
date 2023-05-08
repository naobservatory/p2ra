from pathogen_properties import *

background = """Varicella-zoster virus (VZV) causes both varicella (chickenpox) and herpes zoster (shingles). After a primary infection with VZV,
the virus can remain dormant in the sensory nerve cells of the body,
reactivating later in life and causing herpes zoster, or shingles."""


pathogen_chars = PathogenChars(
    na_type=NAType.DNA,
    enveloped=Enveloped.ENVELOPED,
    taxid=TaxID(10335),
)
